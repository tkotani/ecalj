C ----------------------------------------------------------------
C   spherical Bessel function at x=r(i)*sqrt(e) divided by x**l 
c
Ci Inputs
Ci   r    list of points
Ci   e    energy
Ci   y,h  work vectors of length n each 
Co Outputs
Co   xi   J(r,l)/r**l, according to standard definition
C
c origianl: Feb. 15, 2010, Hiori Kino
      subroutine ropbes(r,e,lmax,y,h,xi,n)
      implicit none
      integer lmax,n
      double precision e,r(n),xi(n,0:lmax),h(n),y(n)

      double precision:: e2,eps=1.0d-10
      integer:: i,l
      e2=sqrt(e)
      do i=1,n
         h(i)= r(i)*e2
         if (h(i)< eps) h(i)=eps 
      enddo
      do l=0,lmax
        call dsjnv(l,n,h,y)
        do i=1,n
        xi(i,l)=y(i)/h(i)**l 
        enddo
      enddo   
 
      end subroutine ropbes

      
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine dsjnv
     >     (N,K,X,
     <      Y)
c                           @(#)dsjnv.f 9.1 97/05/08 14:48:04 
C     Modified from Numpack(DSJ0,DSJ1). Vectorized Version.
C     X should be greater eqoual 0 and less than 3.53d15.
c
c     modified to be able to calculate spherical bessel up to l=7
c                            by H. Sawada   on 12th Jul. 1996
C**********************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER   N,K,I
      REAL*8 X(*),Y(*),T,B,C,D,
     & A00,A01,A02,A03,A04,A05,A06,A07,
     & A10,A11,A12,A13,A14,A15,A16,A17,
     & A20,A21,A22,A23,A24,A25,A26,A27,
     & A30,A31,A32,A33,A34,A35,A36,A37,
     & A40,A41,A42,A43,A44,A45,A46,A47,
     & A50,A51,A52,A53,A54,A55,A56,A57,
     & A60,A61,A62,A63,A64,A65,A66,A67,
     & A70,A71,A72,A73,A74,A75,A76,A77
C N = 0
      DATA A00 /  0.10000000000000000000D+01/
      DATA A01 / -0.16666666666666666667D+00/
      DATA A02 /  0.83333333333333333333D-02/
      DATA A03 / -0.19841269841269841270D-03/
      DATA A04 /  0.27557319223985890653D-05/
      DATA A05 / -0.25052108385441718775D-07/
      DATA A06 /  0.16059043836821614599D-09/
      DATA A07 / -0.76471637318198164759D-12/
C N = 1
      DATA A10 /  0.33333333333333333333D+00/
      DATA A11 / -0.33333333333333333333D-01/
      DATA A12 /  0.11904761904761904762D-02/
      DATA A13 / -0.22045855379188712522D-04/
      DATA A14 /  0.25052108385441718775D-06/
      DATA A15 / -0.19270852604185937519D-08/
      DATA A16 /  0.10706029224547743066D-10/
      DATA A17 / -0.44983316069528332211D-13/
C N = 2
      DATA A20 /  0.66666666666666666667D-01/
      DATA A21 / -0.47619047619047619048D-02/
      DATA A22 /  0.13227513227513227513D-03/
      DATA A23 / -0.20041686708353375020D-05/
      DATA A24 /  0.19270852604185937519D-07/
      DATA A25 / -0.12847235069457291680D-09/
      DATA A26 /  0.62976642497339665096D-12/
      DATA A27 / -0.23675429510278069585D-14/
C N = 3
      DATA A30 /  0.95238095238095238095D-02/
      DATA A31 / -0.52910052910052910053D-03/
      DATA A32 /  0.12025012025012025012D-04/
      DATA A33 / -0.15416682083348750015D-06/
      DATA A34 /  0.12847235069457291680D-08/
      DATA A35 / -0.75571970996807598115D-11/
      DATA A36 /  0.33145601314389297419D-13/
      DATA A37 / -0.11274014052513366469D-15/
C N = 4
      DATA A40 /  0.10582010582010582011D-02/
      DATA A41 / -0.48100048100048100048D-04/
      DATA A42 /  0.92500092500092500093D-06/
      DATA A43 / -0.10277788055565833344D-07/
      DATA A44 /  0.75571970996807598115D-10/
      DATA A45 / -0.39774721577267156903D-12/
      DATA A46 /  0.15783619673518713057D-14/
      DATA A47 / -0.49017452402232028126D-17/
C N = 5
      DATA A50 /  0.96200096200096200096D-04/
      DATA A51 / -0.37000037000037000037D-05/
      DATA A52 /  0.61666728333395000062D-07/
      DATA A53 / -0.60457576797446078492D-09/
      DATA A54 /  0.39774721577267156903D-11/
      DATA A55 / -0.18940343608222455668D-13/
      DATA A56 /  0.68624433363124839376D-16/
      DATA A57 / -0.19606980960892811250D-18/
C N = 6
      DATA A60 /  0.74000074000074000074D-05/
      DATA A61 / -0.24666691333358000025D-06/
      DATA A62 /  0.36274546078467647095D-08/
      DATA A63 / -0.31819777261813725522D-10/
      DATA A64 /  0.18940343608222455668D-12/
      DATA A65 / -0.82349320035749807252D-15/
      DATA A66 /  0.27449773345249935751D-17/
      DATA A67 / -0.72618448003306708335D-20/
C N = 7
      DATA A70 /  0.49333382666716000049D-06/
      DATA A71 / -0.14509818431387058838D-07/
      DATA A72 /  0.19091866357088235313D-09/
      DATA A73 / -0.15152274886577964534D-11/
      DATA A74 /  0.82349320035749807252D-14/
      DATA A75 / -0.32939728014299922901D-16/
      DATA A76 /  0.10166582720462939167D-18/
      DATA A77 / -0.25040844139071278736D-21/

      IF(N.EQ.0) THEN
*POPTION cyclic
          DO 100 I=1,K
            IF(X(I).GT.1.D0) THEN
                Y(I) = SIN(X(I))/X(I)
              ELSE
                T    = X(I)**2
                Y(I) = ((((((A07*T+A06)*T+A05)*T+A04)*T
     &                            +A03)*T+A02)*T+A01)*T+A00
            END IF
  100     CONTINUE
        ELSE IF(N.EQ.1) THEN
*POPTION cyclic

          DO 110 I=1,K
            IF(X(I).GT.1.D0) THEN
                Y(I) = (SIN(X(I))/X(I)-COS(X(I)))/X(I)
              ELSE
                T    = X(I)**2
                Y(I) = (((((((A17*T+A16)*T+A15)*T+A14)*T
     &                             +A13)*T+A12)*T+A11)*T+A10)*X(I)
            END IF
  110     CONTINUE
        ELSE IF(N.EQ.2) THEN
*POPTION cyclic

          DO 120 I=1,K
            IF(X(I).GT.1.D0) THEN
                B    = 1.D0/X(I)
                C    = 3.D0*B
                Y(I) = B*((B*C-1.D0)*SIN(X(I))-C*COS(X(I)))
              ELSE
                T    = X(I)**2
                Y(I) = (((((((A27*T+A26)*T+A25)*T+A24)*T
     &                             +A23)*T+A22)*T+A21)*T+A20)*T
            END IF
  120     CONTINUE
        ELSE IF(N.EQ.3) THEN
*POPTION cyclic

          DO 130 I=1,K
            IF(X(I).GT.1.D0) THEN
                B    = 1.D0/X(I)
                C    = 15.D0*B
                Y(I) = B**2*(  (B*C-6.D0)*SIN(X(I))
     &                        -(C-X(I))  *COS(X(I)) )
              ELSE
                T    = X(I)**2
                Y(I) = (((((((A37*T+A36)*T+A35)*T+A34)*T
     &                             +A33)*T+A32)*T+A31)*T+A30)*T*X(I)
            END IF
  130     CONTINUE
        ELSE IF(N.EQ.4) THEN
*POPTION cyclic

          DO 140 I=1,K
            IF(X(I).GT.1.D0) THEN
                B    = 1.D0/X(I)
                C    = B**2
                D    = 105.D0*C
                Y(I) = B*(  (C*D-45.D0*C+1.D0)*SIN(X(I))
     &                     - B*(D-10.D0)      *COS(X(I)) )
              ELSE
                T    = X(I)**2
                Y(I) = (((((((A47*T+A46)*T+A45)*T+A44)*T
     &                             +A43)*T+A42)*T+A41)*T+A40)*T**2
            END IF
  140     CONTINUE
        ELSE IF(N.EQ.5) THEN
*POPTION cyclic

          DO 150 I=1,K
            IF(X(I).GT.1.D0) THEN
                B    = 1.D0/X(I)
                C    = B**2
                D    = 945.D0*C
                Y(I) = C*(  (C*D-420.D0*C+15.D0) *SIN(X(I))
     &                     - B*(D-105.D0+X(I)**2)*COS(X(I)) )
              ELSE
                T    = X(I)**2
                Y(I) = (((((((A57*T+A56)*T+A55)*T+A54)*T
     &                             +A53)*T+A52)*T+A51)*T+A50)*T**2*X(I)
            END IF
  150     CONTINUE
        ELSE IF(N.EQ.6) THEN
*POPTION cyclic

          DO 160 I=1,K
            IF(X(I).GT.1.D0) THEN
                B    = 1.D0/X(I)
                C    = B**2
                D    = 10395.D0*C*C
                Y(I) = B*(  (C*D-4725.D0*C*C+210.D0*C-1.D0)*SIN(X(I))
     &                     - B*(D-1260.D0*C+21.D0)         *COS(X(I)) )
              ELSE
                T    = X(I)**2
                Y(I) = (((((((A67*T+A66)*T+A65)*T+A64)*T
     &                             +A63)*T+A62)*T+A61)*T+A60)*T**3
            END IF
  160     CONTINUE
        ELSE IF(N.EQ.7) THEN
*POPTION cyclic

          DO 170 I=1,K
            IF(X(I).GT.1.D0) THEN
                B    = 1.D0/X(I)
                C    = B**2
                D    = 135135.D0*C*C
                Y(I) = C*(  (C*D-62370.D0*C*C+3150.D0*C-28.D0)
     &                      *SIN(X(I))
     &                     - B*(D-17325.D0*C+378.D0-X(I)**2)
     &                      *COS(X(I)) )
              ELSE
                T    = X(I)**2
                Y(I) = (((((((A77*T+A76)*T+A75)*T+A74)*T
     &                             +A73)*T+A72)*T+A71)*T+A70)*T**3*X(I)
            END IF
  170     CONTINUE
      END IF
      RETURN
      END
