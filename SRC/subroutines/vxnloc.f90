subroutine pbevxc(nr,nsp,nrx,rofi,rp,exc,vxc,mode,units, &
     grh,ggrh,agrh,grgagr)
  !- Interface for Perdew-Burke-Ernzerhof generalized gradient approximation
  !  10.07.96  dln
  ! ---------------------------------------------------------------
  !i Inputs:
  !i   nr   - No. of points in the radial mesh
  !i   nsp  - 1 spin-restricted, 2- spin-pol.
  !i   nrx  - Max. dimension of arrays (> or = nr)
  !i   rofi - radial points
  !i   rp   - rho
  !i   mode - 0 : LSDA
  !i          1 : PW91
  !i          2 : PBE
  !i         -1 : PW91 correction only
  !i         -2 : PBE  correction only
  !i         10 : LSDA, same as 0, but no final extrapolation
  !i   units - Ry : output in Rydbergs
  !i           not Ry: Hartree
  !o Outputs:
  !o   exc(nrx) -  exchange-corr. energy
  !o   vxc(nrx) -  exchange-corr. potential
  !r Remarks:
  !r   Gradient arrays are NOT used if no GGA is required.
  !r   This saves memory
  !r   working arrays arrays:
  !r   grh -
  !r   ggrh -
  !r   agrh -
  !r   grgagrh -
  !r   References:
  ! ----------------------------------------------------------------
  !     implicit none

  ! - Input variables
  integer :: nr,nsp,nrx,mode,isp
  double precision :: rp(nrx,nsp), &
       exc(nrx),vxc(nrx,nsp),rofi(nr), &
       Xagrh1,Xagrh2,Xgrgagr1,Xgrgagr2,Xggrh1,Xggrh2,XXagrh, &
       XXgrgagr

  ! - Work arrays and internal variables
  double precision :: grh(nrx,2),ggrh(nrx,2), &
       agrh(nrx,4),grgagr(nrx,3)
  double precision :: exlsd, vxuplsd, vxdnlsd, eclsd, &
       vcuplsd, vcdnlsd, expw91, vxuppw91, vxdnpw91, &
       ecpw91, vcuppw91, vcdnpw91, expbe, vxuppbe, &
       vxdnpbe, ecpbe, vcuppbe, vcdnpbe
  integer :: ir,i,iprint
  character*2, units

  if(mode == 0 .OR. mode == 10) GOTO 9999 !! Do not calculate gradients

  call pshpr(iprint()-20)

  ! --- grad(rho), laplacian rho ---
  call radgrx(nr,nrx,nsp,rofi,rp,grh)
  call radgrx(nr,nrx,nsp,rofi,grh,ggrh)
  do  20  i  = 1, nsp
     do  24  ir = 2, nr
        ggrh(ir,i) = ggrh(ir,i) + 2d0*grh(ir,i)/rofi(ir)
24   enddo
     ggrh(1,i) =(rofi(3)*ggrh(2,i)-rofi(2)*ggrh(3,i))/(rofi(3)-rofi(2))

     ! --- grad rho . grad abs grad rho ---
     do  26  ir = 1, nr
        agrh(ir,i) = dabs(grh(ir,i))
26   enddo
     call radgrx(nr,nrx,1,rofi,agrh(1,i),grgagr(1,i))
     do  28  ir = 1, nr
        grgagr(ir,i) = grh(ir,i)*grgagr(ir,i)
28   enddo
20 enddo

  ! --- Extra terms g(n), g(n+).g(n-), g(n).g(abs(g(n))) if spin pol ---
  if (nsp == 2) then
     do  32  ir = 1, nr
        agrh(ir,3) = dabs(grh(ir,1)+grh(ir,2))
32   enddo
     call radgrx(nr,nrx,1,rofi,agrh(1,3),grgagr(1,3))
     do  34  ir = 1, nr
        grgagr(ir,3) = (grh(ir,1)+grh(ir,2))*grgagr(ir,3)
34   enddo
     do  36  ir = 1, nr
        agrh(ir,4) = grh(ir,1)*grh(ir,2)
36   enddo
  endif

  call poppr

  ! --- Call PBE subroutine for each radial point ---
9999 continue
  do  40  ir=1,nr
     if (mode /= 0 .AND. mode /= 10) then
        Xagrh1  =agrh  (ir,1)
        Xgrgagr1=grgagr(ir,1)
        Xggrh1  =ggrh  (ir,1)
        Xagrh2  =agrh  (ir,nsp)
        Xgrgagr2=grgagr(ir,nsp)
        Xggrh2  =ggrh  (ir,nsp)
        XXagrh  =agrh  (ir,2*nsp-1)
        XXgrgagr=grgagr(ir,2*nsp-1)
     else
        Xagrh1  =0d0
        Xgrgagr1=0d0
        Xggrh1  =0d0
        Xagrh2  =0d0
        Xgrgagr2=0d0
        Xggrh2  =0d0
        XXagrh  =0d0
        XXgrgagr=0d0
     endif
     call easypbe(rp    (ir,1),   ! &  rho up
     Xagrh1,         ! &  grad (rho_up)
     Xgrgagr1,       ! &  (grad up).(grad |grad up|)
     Xggrh1,         ! &  Laplacian of up
     rp    (ir,nsp), ! &  rho dwn
     Xagrh2,         ! &  grad (rho_dwn)
     Xgrgagr2,       ! & (grad dwn).(grad |grad dwn|)
     Xggrh2,         ! &  Laplacian of dwn
     XXagrh,         ! &  |grad rho|
     XXgrgagr,       ! &  (grad rho).(grad |grad rho|)
1    ,                ! &  do correlations
1    ,                ! &  do potential
0    ,                ! &  Never do Becke exchange
     exlsd,  ! &  LSD exchange energy density
     vxuplsd,! &  up LSD exchange potential
     vxdnlsd,! &  down LSD exchange potential
     eclsd, ! &  LSD exchange-correlation energy density
     vcuplsd,! &  up LSD exchange-correlation potential
     vcdnlsd,! &  down LSD exchange-correlation potential
     expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, ! & PW91
     expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)       !PBE

     !-  LSDA part -- if mode >= 0 !!! Remember output in Hartree!!
     if(mode == 0 .OR. mode == 10) then
        exc(ir)=exlsd+eclsd
        vxc(ir,1)=vXuplsd+vCuplsd
        if(nsp == 2) vxc(ir,2)=vXdnlsd+vCdnlsd
     endif

     !- PW91 part
     if(mode == 1) then                ! PW91 GGA
        exc(ir)=expw91+ecpw91
        vxc(ir,1)=vxuppw91+vcuppw91
        if(nsp == 2) vxc(ir,2)=vxdnpw91+vcdnpw91
     endif

     !- PBE part
     if (mode == 2) then          !PBE GGA
        exc(ir)=expbe+ecpbe
        vxc(ir,1)=vxuppbe+vcuppbe
        if(nsp == 2) vxc(ir,2)=vxdnpbe+vcdnpbe
     endif

     !- PW correction
     if (mode == -1) then
        exc(ir)=expw91+ecpw91-exlsd-eclsd
        vxc(ir,1)=vxuppw91+vcuppw91-vXuplsd-vCuplsd
        if(nsp == 2) vxc(ir,2)=vxdnpw91+vcdnpw91-vXdnlsd-vCdnlsd
     endif

     !- PBE correction
     if (mode == -2) then
        exc(ir)=expbe+ecpbe-exlsd-eclsd
        vxc(ir,1)=vxuppbe+vcuppbe-vXuplsd-vCuplsd
        if(nsp == 2) vxc(ir,2)=vxdnpbe+vcdnpbe-vXdnlsd-vCdnlsd
     endif

40 enddo

  if (units == 'Ry') then
     do  50  ir=1,nr
        exc(ir)=exc(ir)*2
        do  52  isp=1,nsp
           vxc(ir,isp)=vxc(ir,isp)*2
52      enddo
50   enddo
  endif

  if (mode /= 10) then
     do  60  i = 1, nsp
        vxc(1,i)=(vxc(2,i)*rofi(3)-vxc(3,i)*rofi(2))/(rofi(3)-rofi(2))
60   enddo
  endif
end subroutine pbevxc

