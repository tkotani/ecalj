!> Martix for hilbert transformation, rmat.
! --- this is for KLatexFormula -------------  
! &\text{We use this formula(hilbertmat2):} \int_{a}^{b}\frac{(mw+n)dw}{w-z}=m\Big(z\int_{a}^{b}\frac{dw}{w-z}+(b-a)\Big)+n \int_{a}^{b}\frac{dw}{w-z} \\
! &H(z)=\int_{-w_{max}}^{w_{max}} dw \frac{f(w)}{w-z}\\
! &\approx\sum_i\Big( \int_{L_i}^{C_i}\frac{dw}{w-z} 
! \Big(f(i)\frac{w-C_{i-1}}{C_i-C_{i-1}}+f(i-1)\frac{C_i-w}{C_i-C_{i-1}}\Big)
! +\int_{C_i}^{R_i}\frac{dw}{w-z} 
! \Big(f(i)\frac{C_{i+1}-w}{C_{i+1}-C_{i}}+f(i+1)\frac{w-C_i}{C_{i+1}-C_{i}}\Big)\Big) \ \ : f(w) \text{\ linear interpolation}\\
! &\text{ We assume } f(0)=0, f(N+1)=0, f(-N-1)=0. \\
! &C_{-N-1}=L_{-N} \text{(Left end for negative $w$)},C_{1}=0\text{(Left end for positive $w$)}\\
! &C_{-1}=0 \text{\ \ \ \ \ \ \ (Right end for negative $w$)},C_{N+1}=R_{N}, \text{(Right end for positive $w$)},\\
! &\text{We separate integration for positive $w$ and negative $w$}.
!-------------------------------------------
pure subroutine hilbertmat(zz,nwhis, his_L,his_C,his_R, rmat)
  implicit none
  intent(in):: zz,nwhis,his_L,his_C,his_R
  intent(out):: rmat
  !r   zz is real--->  add  + img*epsz
  !r   zz is complex (and Im(zz)>0) : includes all contribution when Im(zz)>eps
  !o   rmat(-nwhis:nwhis) : rmat(0) is not meaningful.
  !i i-th Histgram bin on real axis are given by [his_L, his_R]. center is his_C.
  !r f(zz) = \int_-x(nwhis)^x(nwhis) f(x)/(x-zz)
  !r       = \sum_{i/=0} rmat(i)*f(i)   ,where f(i) is the average value at i-th bin.
  integer :: iw,nwhis
  real(8),parameter:: eps=1d-8, epsz=1d-13
  complex(8),parameter:: img=(0d0,1d0),  imgepsz =img*epsz 
  real(8):: his_L(-nwhis:nwhis),his_C(-nwhis:nwhis),his_R(-nwhis:nwhis),ci,cim,cip,cicim,cipci
  complex(8)::zz, rr_fac(-nwhis:nwhis),rl_fac(-nwhis:nwhis), domega_c,domega_r,domega_l, rmat(-nwhis:nwhis),facm,facn
!  external:: hilbertmat2
  rmat=0d0
  do iw = -nwhis, nwhis !We assume f(0)=f(nhwis)=f(-nwfis)=0
    if(iw==0) cycle
    ci = his_C(iw)
! \int_Li^Ci  left edge to center. Linear interpolation with f(cim) and f(ci)
    if(iw==-nwhis) then; cim= his_L(iw)
    elseif(iw==1)  then; cim= 0d0
    else               ; cim= his_C(iw-1)
    endif
    cicim = ci-cim
    call hilbertmat2(zz+imgepsz, his_L(iw),his_C(iw), facm,facn)
    rmat(iw)                              = rmat(iw)   + ( facm - facn*cim ) /cicim 
    if(iw/=-nwhis .AND. iw/=1) rmat(iw-1) = rmat(iw-1) + (-facm + facn*ci  ) /cicim 
! \int_Ci^Ri center to right edge. Linear interpolation with f(ci) and f(cip)
    if(iw==nwhis)  then; cip= his_R(iw)
    elseif(iw==-1) then; cip= 0d0
    else               ; cip= his_C(iw+1)
    endif
    cipci=cip-ci
    call hilbertmat2(zz+imgepsz, his_C(iw),his_R(iw), facm,facn)
    rmat(iw)                              = rmat(iw)   + ( -facm + facn*cip ) /cipci 
    if(iw/=nwhis .AND. iw/=-1) rmat(iw+1) = rmat(iw+1) + (  facm - facn*ci  ) /cipci 
 enddo
 contains
   pure subroutine hilbertmat2(zz,a,b, facm,facn) !\int_a^b dw (m*w+n)/(z-w) = m*facm+ n*facn
     intent(in)::              zz,a,b
     intent(out)::                     facm,facn
     real(8):: a,b
     real(8),parameter:: eps=1d-8
     complex(8):: zz,facm,facn,domega_b,domega_a
     domega_b = zz - b 
     domega_a = zz - a
     facm =0d0
     facn =0d0
     if(abs(domega_a)>eps.and.abs(domega_b)>eps) then
       facn =  log(domega_b/domega_a)
       facm =  zz*facn + (b-a) ! b-a is needed for linear interpolation for f(i).
       ! facm =  zz*log(domega_b/domega_a) 
     endif
  end subroutine hilbertmat2
  pure complex(8) function safe_log(z)
    intent(in)::z
    complex(8) :: z,zz
    real(8),parameter:: eps=1d-6
    zz=z
    if(abs(z)<eps    ) zz=z/abs(z)*eps
    if(abs(z)>1d0/eps) zz=z/abs(z)/eps
    safe_log = log(zz)
  end function safe_log
end subroutine hilbertmat

! !> Martix for hilbert transformation, rmat.
! pure subroutine hilbertmat(zz,nwhis, his_L,his_C,his_R, rmat)
!   implicit none
!   intent(in):: zz,nwhis,his_L,his_C,his_R
!   intent(out):: rmat
!   !r  zz is real--->  no img*delta function part
!   !r   zz is complex (and Im(zz)>0) : includes all contribution when Im(zz)>eps
!   !o  rmat(-nwhis:nwhis) : rmat(0) is not meaningful.
!   !i i-th Histgram bin on real axis are given by [his_L, his_R]. center is his_C.
!   !r f(zz) = \int_-x(nwhis)^x(nwhis) f(x)/(zz-x)
!   !r       = \sum_{i/=0} rmat(i)*f(i)   ,where f(i) is the average value at i-th bin.
!   integer :: iw,nwhis
!   real(8),parameter:: eps=1d-8, epsz=1d-13
!   complex(8),parameter:: img=(0d0,1d0)
!   real(8):: his_L(-nwhis:nwhis),his_C(-nwhis:nwhis),his_R(-nwhis:nwhis), delta_r,delta_l,ddr,ddl,delta_rr,delta_ll
!   complex(8)::zz,imgepsz, rr_fac(-nwhis:nwhis),rl_fac(-nwhis:nwhis), domega_c,domega_r,domega_l, rmat(-nwhis:nwhis)
!   imgepsz =img*epsz 
!   rmat=0d0
!   do iw = -nwhis, nwhis
!     if(iw==0) cycle
!     domega_r = zz - his_R(iw) + imgepsz
!     domega_c = zz - his_C(iw) + imgepsz
!     domega_l = zz - his_L(iw) + imgepsz
!     ! rr_fac(his_C(is))=\int^{his_R}_{his_C} d omega' /(his_C(is) -omega')! rr_fac(iw)=log( abs((domega_r/domega_c)) )
!     if( abs(domega_c)<eps .OR. abs(domega_r)<eps ) then; rr_fac(iw) = 0d0
!     else;                                                rr_fac(iw) = log( domega_r/domega_c )
!     endif
!     ! rl_fac(his_C(is)) = \int^{his_C}^{his_L} d omega' /(his_C(is) -omega') !  rl_fac(iw) = log( abs((domega_c/domega_l)) )
!     if( abs(domega_c)<eps .OR. abs(domega_l)<eps ) then; rl_fac(iw) = 0d0 !
!     else ;                                               rl_fac(iw) = log( domega_c/domega_l)
!     endif
    
!     domega_c = zz - his_C(iw)
!     if    (iw==nwhis) then;   delta_rr = his_R(iw)   
!     elseif(iw== -1)   then;   delta_rr =   0d0       
!     else;                     delta_rr = his_C(iw+1) 
!     endif
!     if(iw== -nwhis) then;     delta_ll =   his_L(iw)
!     elseif(iw==  1) then;     delta_ll =   0d0
!     else;                     delta_ll =   his_C(iw-1)
!     endif
!     !          ddr = (his_R(iw)-his_C(iw))/delta_r
!     !          ddl = (his_C(iw)-his_L(iw))/delta_l
! !    rmat(iw)  = rmat(iw  ) + rr_fac(iw)*( 1d0-domega_C/delta_r) !+ ddr
! !    if(iw/=nwhis .AND. iw/=-1) rmat(iw+1) = rmat(iw+1) + rr_fac(iw)*domega_C/delta_r     !- ddr
! !    rmat(iw)  = rmat(iw) + rl_fac(iw)*( 1d0+domega_C/delta_l)   !- ddl
! !    if(iw/=-nwhis .AND. iw/=1) rmat(iw-1) = rmat(iw-1) - rl_fac(iw)*domega_C/delta_l     !+ ddl
    
!     rmat(iw)                              = rmat(iw)   + rl_fac(iw)*(  zz - delta_ll ) /(his_C(iw)-delta_ll)  !- ddl
!     if(iw/=-nwhis .AND. iw/=1) rmat(iw-1) = rmat(iw-1) + rl_fac(iw)*( -zz + his_C(iw)) /(his_C(iw)-delta_ll)  !+ ddl
!     rmat(iw)                              = rmat(iw)   + rr_fac(iw)*( -zz + delta_rr ) /(delta_rr-his_C(iw))  !+ ddr
!     if(iw/=nwhis .AND. iw/=-1) rmat(iw+1) = rmat(iw+1) + rr_fac(iw)*(  zz - his_C(iw)) /(delta_rr-his_C(iw))  !- ddr
!   enddo
! end subroutine hilbertmat

