module m_tetirr
  public tetirr,ccutup
  private
  contains
  subroutine tetirr(qb,n1,n2,n3,ipq,ntet,idtet)
  use m_lmfinit,only: stdo
  !-  Finds inequivalent tetrahedra and counts them
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   ipq   :is a work array of dimension 6*n1*n2*n3.
  !i         :On input, ipq(i1,i2,i3) points to corresponding irreducible qp
  !i   n1..n3:no. of divisions made along each reciprocal lattice vector
  !i   qb    :vectors of first microcell
  !o Outputs:
  !o   idtet :idtet(0,i) = number of tetrahedra of the i'th kind
  !o         :idtet(1-4,i) points to the 4 irreducible k-points defining
  !o         :the tetrahedron.
  !o   ntet  :number of different tetrahedra
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: n1,n2,n3,ntet,ipq(n1,n2,n3),idtet(0:4,*)
  double precision :: qb(3,3)
  integer:: i , j , i1 , i2 , i3 , j1 , j2 , j3 , k1 , k2 , k3 &
       , itet , mtet , ic , ii , ibtr(3,3) , kcut(3,4,6) , imc(0:1,0:1,0:1) &
       , iq(4) , iprint , i1mach
  integer,allocatable :: iwk_iv(:),iprm(:)
  double precision :: qb1(3,3)
  logical :: ipr
  character outs*80
  ipr = iprint() .ge. 30 .and. 6*n1*n2*n3 .gt. 1000.or. iprint() .ge. 40
  if(iprint()>29) write(stdo,"(' TETIRR: sorting ',i8,' tetrahedra ...')")6*n1*n2*n3
  call ccutup(qb,qb1,ibtr,kcut)
  ntet = 0
  ! --- Start loop over microcells ----
  do  202  i3 = 1, n3
     do  201  i2 = 1, n2
        do  20  i1 = 1, n1
           !   ... Set up identifiers at 8 corners of microcell
           do    k1 = 0, 1
              j1 = mod(i1+k1-1,n1) + 1
              do    k2 = 0, 1
                 j2 = mod(i2+k2-1,n2) + 1
                 do    k3 = 0, 1
                    j3 = mod(i3+k3-1,n3) + 1
                    imc(k1,k2,k3) = ipq(j1,j2,j3)
                 enddo
              enddo
           enddo
           !   --- Start loop over tetrahedra ---
           do  10  itet = 1, 6
              do  2  ic = 1, 4
                 k1 = kcut(1,ic,itet)
                 k2 = kcut(2,ic,itet)
                 k3 = kcut(3,ic,itet)
                 iq(ic) = imc(k1,k2,k3)
2             enddo
              !    ... Order the identifiers
              do    j = 1, 3
                 do    i = 1, 4-j
                    if (iq(i) > iq(i+1)) then
                       ii = iq(i)
                       iq(i) = iq(i+1)
                       iq(i+1) = ii
                    endif
                 enddo
              enddo
              ntet = ntet+1
              idtet(0,ntet) = 1
              do  6  i = 1, 4
                 idtet(i,ntet) = iq(i)
6             enddo
10         enddo
20      enddo
201  enddo
202 enddo
  ! --- Eliminate duplicate tetrahedra ---
  mtet = ntet
  allocate(iprm(mtet),iwk_iv(mtet*5))
  call ivheap(5, mtet,idtet,iprm,1)
  call ivprm (5, mtet , idtet , iwk_iv, iprm, .true. )
  deallocate(iwk_iv,iprm)
  ntet = 1
  do  30  i = 2, mtet
     if (idtet(1,i) == idtet(1,ntet)  .AND. &
          idtet(2,i) == idtet(2,ntet)  .AND. &
          idtet(3,i) == idtet(3,ntet)  .AND. &
          idtet(4,i) == idtet(4,ntet)) then
        idtet(0,ntet) = idtet(0,ntet)+1
     else
        ntet = ntet+1
        idtet(1,ntet) = idtet(1,i)
        idtet(2,ntet) = idtet(2,i)
        idtet(3,ntet) = idtet(3,i)
        idtet(4,ntet) = idtet(4,i)
     endif
30 enddo
  if(ipr) write(stdo,"(i8,' inequivalent tetrahedron=')")ntet
end subroutine tetirr
SUBROUTINE CCUTUP(B0,B,IBTR,KCUT)
  IMPLICIT double precision (A-H,O-Z)
  integer:: ibtr,i,ic,itet,j,lx,ly,lz,iprint,lxx,lyy,KCUT0(3,4,6),KCUT(3,4,6)
  DIMENSION B(3,3),IBTR(3,3),B0(3,3),SHIFT(3),P(3,4)
  DATA KCUT0/ &
       0,0,0, 0,1,0, 1,1,0, 1,1,1,  0,0,0, 1,0,0, 1,1,0, 1,1,1, &
       0,0,0, 1,0,0, 1,0,1, 1,1,1,  0,0,0, 0,1,0, 0,1,1, 1,1,1, &
       0,0,0, 0,0,1, 0,1,1, 1,1,1,  0,0,0, 0,0,1, 1,0,1, 1,1,1 /
  DATA SHIFT/0.D0,0.D0,0.D0/
  ANRM2(X,Y,Z)=X*X*1.00001D0+Y*Y*1.00002D0+Z*Z*1.00003D0 &
       -X*0.000004D0-Y*0.000003D0-Z*0.000002D0
  ! ------ CALL CSHEAR TO GET MOST COMPACT CELL (doesn't work) --------
  CALL CSHEAR(B0,B,IBTR)
  ! ----- CHOSE CUTUP WITH SHORTEST MAX EDGE ---------
  if (iprint() > 100) WRITE(*,*) 'CUTUP : '
  LZ = 0
  LXX = 0
  LYY = 0
  EDGMAX = 1.D20
  EDGMIN = 0.D0
  DO  101  LX = 0, 1
     DO  10  LY = 0, 1
        DO  121  ITET = 1, 6
           DO  12  IC = 1, 4
              CALL MXMYMZ(KCUT0(1,IC,ITET),KCUT(1,IC,ITET),LX,LY,LZ)
12         enddo
121     enddo
        EDMIN = 1D20
        EDMAX = 0D0
        DO  20  ITET = 1, 6
           DO  21  IC = 1, 4
              CALL GTBVEC(KCUT(1,IC,ITET),B,SHIFT,P(1,IC))
21         enddo
           DO  131  I = 1, 3
              DO  13  J = I+1, 4
                 XX = ANRM2(P(1,I)-P(1,J),P(2,I)-P(2,J),P(3,I)-P(3,J))
                 EDMAX = DMAX1(EDMAX,XX)
                 EDMIN = DMIN1(EDMIN,XX)
13            enddo
131        enddo
20      enddo
        if (iprint() > 100) &
             WRITE(*,706) LX,LY,DSQRT(EDMIN),DSQRT(EDMAX)
706     FORMAT(' LX,LY=',2I5,'   EDMIN=',F10.5,'   EDMAX=',F10.5)
        IF (EDMAX < EDGMAX) THEN
           LXX = LX
           LYY = LY
           EDGMAX = EDMAX
           EDGMIN = EDMIN
        ENDIF
10   enddo
101 enddo
  DO  221  ITET = 1, 6
     DO  22  IC = 1, 4
        CALL MXMYMZ(KCUT0(1,IC,ITET),KCUT(1,IC,ITET),LXX,LYY,LZ)
22   enddo
221 enddo
  if (iprint() > 100) &
       WRITE(*,783) LXX,LYY,DSQRT(EDGMIN),DSQRT(EDGMAX)
783 FORMAT(' LXX=',I1,'   LYY=',I1,'   EDGMIN=',F10.5, '   EDGMAX=',F10.5)
END SUBROUTINE CCUTUP
SUBROUTINE CSHEAR(B0,B,IBTR)
  !- Tries to make microcell more compact by shearing.
  ! ----------------------------------------------------------------
  !i Inputs
  !i
  !o Outputs
  !o
  !r Remarks
  !r   Ibtr gives the transformation from bo to b.
  ! ----------------------------------------------------------------
  !      implicit none
  ! Passed parameters
  ! Local parameters
  IMPLICIT double precision (A-H,O-Z)
  integer:: ibtr,i,j
  DIMENSION B(3,3),IBTR(3,3),B0(3,3)
  !     logical print
  DO  8  I = 1, 3
     DO  7  J = 1, 3
        B(J,I) = B0(J,I)
        IBTR(J,I) = 0
7    enddo
     IBTR(I,I) = 1
8 enddo
END SUBROUTINE CSHEAR

SUBROUTINE MXMYMZ(KIN,K,LX,LY,LZ)
  !- Do mirrors in x,y,z if lx,ly,lz=1, respectively
  ! ----------------------------------------------------------------
  !i Inputs
  !i
  !o Outputs
  !o
  !r Remarks
  !r
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: KIN(3),K(3),lx,ly,lz
  K(1) = KIN(1)
  K(2) = KIN(2)
  K(3) = KIN(3)
  IF (LX == 1) K(1) = 1-K(1)
  IF (LY == 1) K(2) = 1-K(2)
  IF (LZ == 1) K(3) = 1-K(3)
  RETURN
END SUBROUTINE MXMYMZ

SUBROUTINE GTBVEC(K,B,SHIFT,V)
  implicit none
  integer :: K(3)
  double precision :: V(3),B(3,3),SHIFT(3)
  V(1) = SHIFT(1) + K(1)*B(1,1) + K(2)*B(1,2) + K(3)*B(1,3)
  V(2) = SHIFT(2) + K(1)*B(2,1) + K(2)*B(2,2) + K(3)*B(2,3)
  V(3) = SHIFT(3) + K(1)*B(3,1) + K(2)*B(3,2) + K(3)*B(3,3)
END SUBROUTINE GTBVEC
end module m_tetirr
