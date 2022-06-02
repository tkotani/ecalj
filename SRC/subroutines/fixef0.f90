      subroutine fixef0(zval,nsp,nspc,nevx,ndev,evl,dosw,ef0)
      use m_lgunit,only:stdo
C- Corrects estimate for Fermi level
C ----------------------------------------------------------------------
Ci Inputs
Ci   zval  :valence charge
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nevx  :max number of eigenvalues calculated
Ci   ndev  :leading dimension of evl
Ci   evl   :eigenvalues
Ci   dosw  :dos window
Cio Inputs/Outputs
Cio  ef0   :on input, estimate for Fermi energy
Cio        :on output, revised estimate, if ef0 outside bounds
Cio  dosw  :on input dos window
Cio        :on output, revised if ebot<dosw(1) or dosw(2)<ef0
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,nspc,ndev,nevx
      double precision zval,ef0,evl(ndev*2),dosw(2)!takao evl(ndev)->evl(ndev*2)
C ... Local parameters
      integer:: i , i1 , ipr , nbpw , i1mach
      integer ,allocatable :: bmap_iv(:)
      real(8) ,allocatable :: wk_rv(:)

      double precision w2,xx,doso(2)
C ... Heap

      call getpr(ipr)
c      stdo = lgunit(1)
      if (nsp .eq. 2) then
        nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
        allocate(bmap_iv(abs(-(nevx*nsp/nbpw+1))))
        if (-(nevx*nsp/nbpw+1)<0) bmap_iv(:)=0

        allocate(wk_rv(nevx*nsp))

        call ebcpl ( 0 , ndev , nevx , nsp , nspc , 1 , nbpw , bmap_iv 
     .  , wk_rv , evl )

      endif
      i = max(1,int(zval)/(3-nsp))
      if (ef0 .lt. evl(i)) then
        i1 = (zval + 0.001d0)/(3-nsp)
        w2 = zval/(3-nsp)-i1
        xx = (1-w2)*evl(i1)+w2*evl(i1+1)
        if (ipr.ge.10) write(stdo,"(' Est Ef = ',f15.8,' < evl(',f15.8,i0,')=',f15.8,
     .  ' ... using qval=',f15.8,' revise to ',f15.8)") ef0,i,evl(i),zval,xx
        ef0 = xx
      endif

      if (nsp .eq. 2) then
        call ebcpl ( 1 , ndev , nevx , nsp , nspc , 1 , nbpw , bmap_iv 
     .  , wk_rv , evl )

        if (allocated(wk_rv)) deallocate(wk_rv)
        if (allocated(bmap_iv)) deallocate(bmap_iv)

      endif

      if (dosw(1) .gt. evl(1) .or. dosw(2) .lt. ef0) then
        doso(1) = dosw(1)
        doso(2) = dosw(2)
        dosw(1) = evl(1) - 0.5d0
        dosw(2) = ef0  + 0.5d0
        if (ipr.ge.10) write(stdo,"(' DOS window (',f15.7,x,f15.7,')',
     .  ' reset to (',f15.7,x,f15.7,')')") doso(1),doso(2),dosw(1),dosw(2)
      endif

      end









