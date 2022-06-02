subroutine pkl2ro(mode,ib,rsm,kmax,nr,nlml,nsp,rofi,rwgt,k0,nlm0, &
     fklc,fklr,rho1,rho2,qmx)
  !- Put PkL or GkL expansion of a function on a radial mesh for one site
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :a compound of digits :
  !i         :1s digit
  !i         :  0 add P_kL expansion of rho to rho1 (and possibly rho2)
  !i         :    Here, fkl are coffs to P_kL expansion
  !i         :  1 add G_kL expansion of rho to rho1 (and possibly rho2)
  !i         :    Here, fkl are coffs to G_kL expansion
  !i         :10s digit
  !i         :  0 add expansion to rho1 only; rho2 is not touched
  !i         :  1 add expansion to both rho1 and rho2
  !i         :100s digit
  !i         :  0 coefficients fkl are real (uses fklr)
  !i         :  1 coefficients fkl are complex (uses fklc)
  !i         :1000s digit
  !i         :  0 do not initialize rho1,rho2 before adding expansion
  !i         :  1 initialize rho1,rho2 before adding expansion
  !i         :10000s digit for spin polarized case
  !i         :  1 zero out charge, rho1+ + rho1- (and rho2+ + rho2-)
  !i         :  2 zero out spin, rho1+ - rho1- (and rho2+ + rho2-)
  !i   ib    :site index (used in addressing fkl expansion)
  !i   rsm   :smoothing radius for P_kL (or G_kL) expansion
  !i   kmax  :k-cutoff for P_kL (or G_kL) expansion
  !i   nr    :number of radial mesh points
  !i   nlml  :L-cutoff for rho1,rho2
  !i   nsp   :2 if fkl is spin-pol; 1 if not
  !i   rofi  :radial mesh points
  !i   rwgt  :radial mesh weights
  !i   k0    :leading dimension of fkl
  !i   nlm0  :second dimension of fkl
  !i   fklc  :complex coefficients to P_kL (or G_kL) expansion for this site
  !i         :(only one of fklc or fklr is used; see mode)
  !i   fklr  :real coefficients to P_kL (or G_kL) expansion for this site
  !i         :(only one of fklc or fklr is used; see mode)
  !o Outputs
  !o  rho1  :fkL P_kL (or G_kL) added to local true density rho1
  !o  rho2  :fkL P_kL (or G_kL) added to local smoothed density rho2
  !i   qmx   :charge in rho1 after fkl PkL is added
  !r Remarks
  !u Updates
  !u   21 Nov 01 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,ib,k0,kmax,nlm0,nr,nlml,nsp
  double precision :: qmx,rsm
  double precision :: rofi(nr),rwgt(nr)
  double precision :: rho1(nr,nlml,nsp),rho2(nr,nlml,nsp)
  double complex   fklc(0:k0,nlm0,nsp,ib)
  double precision :: fklr(0:k0,nlm0,nsp,ib)
  ! ... Local parameters
  integer :: lmx,i,ilm,isp,k,l,ll,lmxl,mode0,mode1,mode2,mode3,mode4,np
  double precision :: add,pi,r,rl,srfpi,sum1
  parameter(lmx=10)
  double precision :: fkl(0:kmax,nlm0,2),pkl(0:kmax,0:lmx)

  pi = 4d0*datan(1d0)
  srfpi = dsqrt(4*pi)
  lmxl = ll(nlml)
  mode0 = mod(mode,10)
  mode1 = mod(mode/10,10)
  mode2 = mod(mode/100,10)
  mode3 = mod(mode/1000,10)
  mode4 = mod(mode/10000,10)

  if (lmxl > lmx) call rxi('pklr2o: increase lmx, need',lmxl)
  if (nlml > nlm0) call rxi('pklr2o: increase nlm0, need',nlml)

  do  isp = 1, nsp
     do  ilm = 1, nlml
        do  k = 0, kmax
           if (mode2 == 0) then
              fkl(k,ilm,isp) = fklr(k,ilm,isp,ib)
           else
              fkl(k,ilm,isp) = dble(fklc(k,ilm,isp,ib))
           endif
        enddo
     enddo
  enddo

  if (mode3 == 1) call dpzero(rho1,nr*nlml*nsp)
  if (mode3 /= 0 .AND. mode1 /= 0) call dpzero(rho2,nr*nlml*nsp)

  do  i = 2, nr
     r = rofi(i)
     if (mode0 == 0) then
        call radpkl(r,rsm,kmax,lmxl,kmax,pkl)
     else
        call radgkl(r,rsm,kmax,lmxl,kmax,pkl)
     endif
     do  ilm = 1, nlml
        l = ll(ilm)
        rl = r**l
        do  isp = 1, nsp
           if (mode1 /= 0) then
              do  k = 0, kmax
                 add = fkl(k,ilm,isp)*pkl(k,l)*r*r*rl
                 rho1(i,ilm,isp) = rho1(i,ilm,isp) + add
                 rho2(i,ilm,isp) = rho2(i,ilm,isp) + add
              enddo
           else
              do  k = 0, kmax
                 add = fkl(k,ilm,isp)*pkl(k,l)*r*r*rl
                 rho1(i,ilm,isp) = rho1(i,ilm,isp) + add
              enddo
           endif
        enddo
     enddo
  enddo
  sum1 = 0d0
  do  isp = 1, nsp
     do  i = 1, nr
        sum1 = sum1 + rwgt(i)*rho1(i,1,isp)
     enddo
  enddo
  qmx = sum1*srfpi

  ! zero out spin or charge
  if (mode4 == 0 .OR. nsp /= 2) return
  i = 20                   ! No core
  if (mode1 == 0) i = 30 ! No rho2
  call splrho(i,nsp,nr,nlml,rho1,rho2,sum1)
  np = nr*nlml
  if (mode4 == 1) then   ! Zero density
     call dpzero(rho1(1,1,1),np)
     if (mode1 /= 0) then ! Including rho2
        call dpzero(rho2(1,1,1),np)
     endif
  endif
  if (mode4 == 2) then   ! Zero spin
     call dpzero(rho1(1,1,nsp),np)
     if (mode1 /= 0) then ! Including rho2
        call dpzero(rho2(1,1,nsp),np)
     endif
  endif
  call splrho(i+1,nsp,nr,nlml,rho1,rho2,sum1)

end subroutine pkl2ro

