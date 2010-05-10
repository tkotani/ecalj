Cdelw1 var ['obmap', '1', 'obmap', '-(nevx*nsp/nbpw+1)', 'integer', 'iv_w_', 'fixef0']
Cdelw1 var ['owk', '1', 'owk', 'nevx*nsp', 'real(8)', 'rv_w_', 'fixef0']
      subroutine fixef0(zval,nsp,nspc,nevx,ndev,evl,dosw,ef0)
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
      double precision zval,ef0,evl(ndev),dosw(2)
C ... Local parameters
Cdelw1       integer i,i1,stdo,lgunit,ipr,nbpw,i1mach,obmap,owk
       integer:: i , i1 , stdo , lgunit , ipr , nbpw , i1mach 
       integer ,allocatable :: iv_w_obmap(:)
       real(8) ,allocatable :: rv_w_owk(:)

      double precision w2,xx,doso(2)
C ... Heap
      integer w(1)
      common /w/ w

      call getpr(ipr)
      stdo = lgunit(1)
      if (nsp .eq. 2) then
        nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
Cdelw1         call defi(obmap,-(nevx*nsp/nbpw+1))
         allocate(iv_w_obmap(abs(-(nevx*nsp/nbpw+1))))
         if (-(nevx*nsp/nbpw+1)<0) iv_w_obmap(:)=0

Cdelw1         call defdr(owk, nevx*nsp)
         allocate(rv_w_owk(nevx*nsp))
         if (nevx*nsp<0) rv_w_owk(:)=0.0d0

Cdelw1         call ebcpl(0,ndev,nevx,nsp,nspc,1,nbpw,w(obmap),w(owk),evl)
         call ebcpl ( 0 , ndev , nevx , nsp , nspc , 1 , nbpw , iv_w_obmap 
     .   , rv_w_owk , evl ) 

      endif
      i = max(1,int(zval)/(3-nsp))
      if (ef0 .lt. evl(i)) then
        i1 = (zval + 0.001d0)/(3-nsp)
        w2 = zval/(3-nsp)-i1
        xx = (1-w2)*evl(i1)+w2*evl(i1+1)
        if (ipr.ge.10) call awrit5(' Est Ef = %,3;3d < evl(%i)='//
     .    '%,3;3d ... using qval=%,1;1d, revise to %,4;4d',
     .      ' ',80,stdo,ef0,i,evl(i),zval,xx)
          ef0 = xx
      endif

      if (nsp .eq. 2) then
Cdelw1         call ebcpl(1,ndev,nevx,nsp,nspc,1,nbpw,w(obmap),w(owk),evl)
         call ebcpl ( 1 , ndev , nevx , nsp , nspc , 1 , nbpw , iv_w_obmap 
     .   , rv_w_owk , evl ) 

Cdelw1 rlse name= obmap old_list= obmap owk 
Cdelw1 rlse name= obmap new_list= (None)
Cdelw1         call rlse(obmap)
         if (allocated(rv_w_owk)) deallocate(rv_w_owk)
         if (allocated(iv_w_obmap)) deallocate(iv_w_obmap)

      endif

      if (dosw(1) .gt. evl(1) .or. dosw(2) .lt. ef0) then
        doso(1) = dosw(1)
        doso(2) = dosw(2)
        dosw(1) = evl(1) - 0.5d0
        dosw(2) = ef0  + 0.5d0
        if (ipr.ge.10) call awrit4(' (warning) DOS window (%;3d,%;3d)'//
     .    ' reset to (%,4;4d,%,4;4d)',' ',80,stdo,doso(1),doso(2),
     .    dosw(1),dosw(2))
      endif

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end









