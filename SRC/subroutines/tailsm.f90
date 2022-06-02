subroutine tailsm(lrhot,nr,nrmt,nsp,a,b,rmt,rsm,nxi0,nxi,exi,rofi, &
     rho,rhot,hfc,hfct)
  use m_lgunit,only:stdo
  !- Fit tails of rho to smoothed Hankel functions
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lrhot :0 make fit for rho only; 1 make fit for rho and rhot
  !i   nr    :number of radial mesh points
  !i   nrmt  :number of points between 0..rmt
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   rmt   :muffin-tin radius, in a.u.
  !i   rsm   :smoothing radius for smoothed hankel fit
  !i   nxi0  :leading dimension of hfc, nfct
  !i   nxi   :number of energies to include in fit
  !i   exi   :energies to include in fit
  !i   rofi  :radial mesh points
  !i   rho   :spherical valence charge density times 4*pi*r*r
  !i   rhot  :total charge density (used if lrhot=1)
  !o Outputs
  !o   hfc   :coefficients to h.f. fits to rho
  !o   hfct  :coefficients to h.f. fits to rhot (if lrhot=1)
  !l Local variables
  !l   rsq   :work array holding rofi**2
  !r Remarks
  !r   A fit is make to tails of the valence charge density for r>rmt
  !r   using smoothed hankel functions of smoothing radius rsm.
  !r   Fit is constrained to agree with integrated charge.
  !u Updates
  !u   19 Apr 02 Move rsq to a local array.  Altered argument list.
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lrhot,nsp,nr,nrmt,nxi0,nxi
  double precision :: a,b,rmt,rsm
  double precision :: hfc(nxi0,nsp),exi(nxi),rho(*),rofi(nr)
  double precision :: hfct(nxi0,nsp),rhot(*)
  ! ... Local parameters
  logical :: lzero
  integer:: lx0(20) , isp , ie , ir , ipr
  !      real(8) ,allocatable :: wk_rv(:)
  real(8) ,allocatable :: idx_rv(:)
  real(8) ,allocatable :: xi_rv(:)

  double precision :: sr,x0(0:2),xi(0:2),fpi,qout,qcst(20),qcst0(20), &
       err,rsq(nr)
  ! ... Heap
  print *,'tailsm: init'
  fpi = 16*datan(1d0)
  call getpr(ipr)
  !      stdo = lgunit(1)

  ! --- Tabulate smoothed Hankels on a radial mesh ---
  do  10  ir = 1, nr
     rsq(ir) = rofi(ir)**2
10 enddo
  do  12  ie = 1, nxi
     lx0(ie) = 0
12 enddo
  allocate(xi_rv(nr*nxi))

  allocate(idx_rv(nr*2))

  !      allocate(wk_rv(nr*6))

  !     Patch to properly handle case when rsm=0
  lzero = rsq(1) .eq. 0 .and. rsm .eq. 0
  if (lzero) rsq(1) = rsq(2)/100
  !      print *,'tailsm:xxx1'
  !      call hansr ( rsm , 0 , 0 , nxi , lx0 , exi , rsq , nr , nr ,
  !     .idx_rv , wk_rv , 0 , xi_rv )
  call hansr ( rsm , 0 , 0 , nxi , lx0 , exi , rsq , nr , nr , &
       idx_rv, 0 , xi_rv )

  !      print *,'tailsm:xxx2'
  if (lzero) rsq(1) = 0
  !      if (allocated(wk_rv)) deallocate(wk_rv)
  if (allocated(idx_rv)) deallocate(idx_rv)


  ! --- Fit smoothed hankels to rho ---
  do  20  isp = 1, nsp

     if (isp == 1 .AND. ipr >= 30) write(stdo,333) nxi,rmt,rsm
333  format(/' tailsm: fit tails to',i2,' smoothed hankels, rmt=', &
          f8.5,', rsm=',f8.5)
     if (isp == 2 .AND. ipr >= 20) &
          write(stdo,'(/'' tailsm: spin 2 ...'')')

     !   ... Fitting constraints for smoothed Hankels
     do  14  ie = 1, nxi
        !          print *,'eee exi=',ie,exi(ie)
        if (rsm < 1d-9) then
           sr = dsqrt(-exi(ie))*rmt
           qcst(ie) = -dsqrt(fpi)*(sr+1)*dexp(-sr)/exi(ie)
        else
           call hansmr(rmt,0d0,1/rsm,x0,1)
           call hansmr(rmt,exi(ie),1/rsm,xi,1)
           qcst(ie) = dsqrt(fpi)/exi(ie)*(-dexp(rsm**2/4*exi(ie)) &
                - rmt**3*(xi(1)-dexp(rsm**2/4*exi(ie))*x0(1)))
        endif
        qcst0(ie) = qcst(ie)
14   enddo

     !   ... Fit for this spin
     call hnsmft ( rofi , rho ( 1 + ( isp - 1 ) * nr ) , nr , qout &
          , a , b , nrmt , exi , qcst , xi_rv , hfc ( 1 , isp ) , nxi &
          , err )


     !        if (isp .eq. 1 .and. ipr .ge. 20 .and. ipr .lt. 30)
     !     .  call awrit3('%N tailsm:  fit tails to %i functions with'
     !     .  //' rsm=%;7g.  rms error=%;7F',' ',80,stdo,nxi,rsm,err)
     if (isp==1 .AND. ipr>=20)  write(stdo,"(' tailsm:  fit tails to ',i8,' functions with' &
          //' rsm=',d13.5,' rms error=',d13.5)") nxi,rsm,err

     !   ... Fit a second time for the full density in rhot
     if (lrhot /= 0) then
        do  15  ie=1,nxi
           qcst(ie)=qcst0(ie)
15      enddo
        call hnsmft ( rofi , rhot ( 1 + ( isp - 1 ) * nr ) , nr , qout &
             , a , b , nrmt , exi , qcst , xi_rv , hfct ( 1 , isp ) , nxi &
             , err )

     endif

     !C   ... Evaluate integral of smoothed charge density (printout)
     !        qsm = 0
     !        qsmr = 0
     !        qsmt = 0
     !        qsmrt = 0
     !        rmt0 = 0
     !        do  16  ie = 1, nxi
     !          if (rsm .lt. 1d-9) then
     !            sr = dsqrt(-exi(ie))*rmt0
     !            xx = -dsqrt(fpi)*(sr+1)*dexp(-sr)/exi(ie)
     !          else
     !            call hansmr(rmt0,0d0,1/rsm,x0,1)
     !            call hansmr(rmt0,exi(ie),1/rsm,xi,1)
     !            xx = dsqrt(fpi)/exi(ie)*(-dexp(rsm**2/4*exi(ie))
     !     .        - rmt0**3*(xi(1)-dexp(rsm**2/4*exi(ie))*x0(1)))
     !          endif
     !          qsmr = qsmr + hfc(ie,isp)*qcst0(ie)
     !          qsm  = qsm  + hfc(ie,isp)*xx
     !          qsmrt = qsmrt + hfct(ie,isp)*qcst0(ie)
     !          qsmt  = qsmt  + hfct(ie,isp)*xx
     !c|          write (stdo,455) ie,xx,qcst0(ie)
     !c|  455     format('  ie=',i4,'   qall',f14.8,'   qout',f14.8)

     !   16   continue
     !        if (ipr .ge. 30) print 345, qsm-qsmr,qsm,qsmt-qsmrt,qsmt
     !  345   format (' valence:   qsm(r<rmt)=',f10.6,
     !     .     '  qsm(all space)=',f10.6/
     !     .     ' val+core:  qsm(r<rmt)=',f10.6,'  qsm(all space)=',f10.6)



20 enddo
  if (allocated(xi_rv)) deallocate(xi_rv)


  print *,'tailsm: end'
end subroutine tailsm

