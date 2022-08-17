subroutine tailsm(lrhot,nr,nrmt,nsp,a,b,rmt,rsm,nxi0,nxi,exi,rofi, rho,rhot,hfc,hfct)
  !- Fit tails of rho to smoothed Hankel functions
  use m_lgunit,only:stdo
  use m_hansr,only :hansr,hansmr
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
  implicit none
  integer :: lrhot,nsp,nr,nrmt,nxi0,nxi
  double precision :: a,b,rmt,rsm
  double precision :: hfc(nxi0,nsp),exi(nxi),rho(*),rofi(nr),hfct(nxi0,nsp),rhot(*)
  logical :: lzero
  integer:: lx0(20) , isp , ie , ir , ipr
  integer,allocatable :: idx_rv(:)
  real(8) ,allocatable :: xi_rv(:)
  double precision :: sr,x0(0:2),xi(0:2),fpi,qout,qcst(20),qcst0(20), err,rsq(nr)
  print *,'tailsm: init'
  fpi = 16*datan(1d0)
  call getpr(ipr)
  ! --- Tabulate smoothed Hankels on a radial mesh ---
  rsq = rofi**2
  lx0=0d0
  allocate(xi_rv(nr*nxi))
  allocate(idx_rv(nr*2))
  !     Patch to properly handle case when rsm=0
  lzero = rsq(1) .eq. 0 .and. rsm .eq. 0
  if (lzero) rsq(1) = rsq(2)/100
  call hansr ( rsm , 0 , 0 , nxi , lx0 , exi , rsq , nr , nr , idx_rv, 0 , xi_rv )
  if (lzero) rsq(1) = 0
  deallocate(idx_rv)
  ! --- Fit smoothed hankels to rho ---
  isploop: do  20  isp = 1, nsp
     if (isp == 1 .AND. ipr >= 30) write(stdo,333) nxi,rmt,rsm
333  format(/' tailsm: fit tails to',i2,' smoothed hankels, rmt=', f8.5,', rsm=',f8.5)
     if (isp == 2 .AND. ipr >= 20) write(stdo,'(/'' tailsm: spin 2 ...'')')
     !   ... Fitting constraints for smoothed Hankels
     do  ie = 1, nxi
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
     enddo
     !   ... Fit for this spin
     call hnsmft ( rofi , rho ( 1 + ( isp - 1 ) * nr ) , nr , qout &
          , a , b , nrmt , exi , qcst , xi_rv , hfc ( 1 , isp ) , nxi , err )
     if (isp==1 .AND. ipr>=20)  write(stdo,"(' tailsm:  fit tails to ',i8,' functions with' &
          //' rsm=',d13.5,' rms error=',d13.5)") nxi,rsm,err
     !   ... Fit a second time for the full density in rhot
     if (lrhot /= 0) then
        do  ie=1,nxi
           qcst(ie)=qcst0(ie)
        enddo
        call hnsmft ( rofi , rhot ( 1 + ( isp - 1 ) * nr ) , nr , qout &
             , a , b , nrmt , exi , qcst , xi_rv , hfct ( 1 , isp ) , nxi , err )
     endif
20 enddo isploop
  if (allocated(xi_rv)) deallocate(xi_rv)
end subroutine tailsm
