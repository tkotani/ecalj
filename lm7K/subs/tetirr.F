Cdelw1 var ['oiwk', '1', 'oiwk', 'mtet*5', 'integer', 'iv_w_', 'tetirr']
      subroutine tetirr(qb,n1,n2,n3,ipq,ntet,idtet)
C-  Finds inequivalent tetrahedra and counts them
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ipq   :is a work array of dimension 6*n1*n2*n3.
Ci         :On input, ipq(i1,i2,i3) points to corresponding irreducible qp
Ci         :ipq is DESTROYED on output.
Ci   n1..n3:no. of divisions made along each reciprocal lattice vector
Ci   qb    :vectors of first microcell
Co Outputs:
Co   idtet :idtet(0,i) = number of tetrahedra of the i'th kind
Co         :idtet(1-4,i) points to the 4 irreducible k-points defining
Co         :the tetrahedron.
Co   ntet  :number of different tetrahedra
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer n1,n2,n3,ntet,ipq(n1,n2,n3),idtet(0:4,*)
      double precision qb(3,3)
C Local variables:
Cdelw1       integer i,j,i1,i2,i3,j1,j2,j3,k1,k2,k3,itet,mtet,ic,ii,ibtr(3,3),
Cdelw1      .  kcut(3,4,6),imc(0:1,0:1,0:1),iq(4),iprint,awrite,i1mach,oiwk
       integer:: i , j , i1 , i2 , i3 , j1 , j2 , j3 , k1 , k2 , k3 
     . , itet , mtet , ic , ii , ibtr(3,3) , kcut(3,4,6) , imc(0:1,0:1,0:1) 
     . , iq(4) , iprint , awrite , i1mach 
       integer ,allocatable :: iv_w_oiwk(:)

      double precision qb1(3,3)
      logical ipr
      character outs*80
C Heap:
      integer w(1)
      common /w/ w

      ipr = iprint() .ge. 30 .and. 6*n1*n2*n3 .gt. 1000
     .  .or. iprint() .ge. 40
      i3 = awrite(' TETIRR: sorting %i tetrahedra ...',outs,80,0,
     .  6*n1*n2*n3,0,0,0,0,0,0,0)
      if (ipr) call cwrite(outs,0,i3-1,0)
      call ccutup(qb,qb1,ibtr,kcut)
      ntet = 0
C --- Start loop over microcells ----
      do  20  i3 = 1, n3
        do  20  i2 = 1, n2
        do  20  i1 = 1, n1
C   ... Set up identifiers at 8 corners of microcell
        do  1  k1 = 0, 1
          j1 = mod(i1+k1-1,n1) + 1
          do  1  k2 = 0, 1
          j2 = mod(i2+k2-1,n2) + 1
          do  1  k3 = 0, 1
          j3 = mod(i3+k3-1,n3) + 1
          imc(k1,k2,k3) = ipq(j1,j2,j3)
    1   continue
C   --- Start loop over tetrahedra ---
        do  10  itet = 1, 6
          do  2  ic = 1, 4
            k1 = kcut(1,ic,itet)
            k2 = kcut(2,ic,itet)
            k3 = kcut(3,ic,itet)
            iq(ic) = imc(k1,k2,k3)
    2     continue
C    ... Order the identifiers
          do  3  j = 1, 3
          do  3  i = 1, 4-j
            if (iq(i) .gt. iq(i+1)) then
              ii = iq(i)
              iq(i) = iq(i+1)
              iq(i+1) = ii
            endif
    3     continue
          ntet = ntet+1
          idtet(0,ntet) = 1
          do  6  i = 1, 4
            idtet(i,ntet) = iq(i)
    6     continue
   10   continue
   20 continue

C --- Eliminate duplicate tetrahedra ---
      mtet = ntet
C ... Sort it with fast ivheap unless there is not enough memory
      call defask(i)
      if (i .gt. mtet*5) then
        call ivheap(5,mtet,idtet,ipq,1)
Cdelw1         call defi(oiwk,mtet*5)
         allocate(iv_w_oiwk(mtet*5))
Ctakao_ZeroClear_NotRequiered          if (mtet*5<0) iv_w_oiwk(:)=0

Cdelw1         call ivprm(5,mtet,idtet,w(oiwk),ipq,.true.)
         call ivprm ( 5 , mtet , idtet , iv_w_oiwk , ipq , .true. ) 

Cdelw1 rlse name= oiwk old_list= oiwk 
Cdelw1 rlse name= oiwk new_list= (None)
Cdelw1         call rlse(oiwk)
         if (allocated(iv_w_oiwk)) deallocate(iv_w_oiwk)

*       call yprm('idtet',0,idtet,0,5,5,mtet)
      else
        call ivshel(5,mtet,idtet,ipq,.false.)
*       call yprm('idtet',0,idtet,0,5,5,mtet)
      endif
      ntet = 1
      do  30  i = 2, mtet
C        print 337, i, idtet(1,i), idtet(2,i), idtet(3,i), idtet(4,i)
C  337   format(i4,4i6)
        if (idtet(1,i) .eq. idtet(1,ntet)  .and.
     .      idtet(2,i) .eq. idtet(2,ntet)  .and.
     .      idtet(3,i) .eq. idtet(3,ntet)  .and.
     .      idtet(4,i) .eq. idtet(4,ntet)) then
          idtet(0,ntet) = idtet(0,ntet)+1
        else
          ntet = ntet+1
          idtet(1,ntet) = idtet(1,i)
          idtet(2,ntet) = idtet(2,i)
          idtet(3,ntet) = idtet(3,i)
          idtet(4,ntet) = idtet(4,i)
        endif
   30 continue

      if (ipr) call awrit1(' %i inequivalent ones found',' ',
     .  80,i1mach(2),ntet)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

