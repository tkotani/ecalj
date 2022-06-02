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
