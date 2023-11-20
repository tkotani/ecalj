!> Martix for hilbert transformation, rmat.
subroutine hilbertmat(zz,nwhis, his_L,his_C,his_R, rmat) 
  !r  zz is real--->  no img*delta function part
  !r   zz is complex (and Im(zz)>0) : includes all contribution when Im(zz)>eps
  !o  rmat(-nwhis:nwhis) : rmat(0) is not meaningful.
  !i i-th Histgram bin on real axis are given by [his_L, his_R]. center is his_C.
  !r f(zz) = \int_-x(nwhis)^x(nwhis) f(x)/(zz-x)
  !r       = \sum_{i/=0} rmat(i)*f(i)   ,where f(i) is the average value at i-th bin.
  implicit none
  integer :: iw,nwhis
  real(8):: his_L(-nwhis:nwhis),his_C(-nwhis:nwhis),his_R(-nwhis:nwhis), eps=1d-8, epsz=1d-13,delta_r,delta_l,ddr,ddl
  complex(8)::zz,imgepsz, rr_fac(-nwhis:nwhis),rl_fac(-nwhis:nwhis),img=(0d0,1d0), domega_c,domega_r,domega_l, rmat(-nwhis:nwhis)
  imgepsz =img*epsz
  rmat=0d0
  do iw = -nwhis, nwhis
     if(iw==0) cycle
     domega_r = zz - his_R(iw) + imgepsz
     domega_c = zz - his_C(iw) + imgepsz
     domega_l = zz - his_L(iw) + imgepsz
     ! rr_fac(his_C(is))=\int^{his_R}_{his_C} d omega' /(his_C(is) -omega')! rr_fac(iw)=log( abs((domega_r/domega_c)) )
     if( abs(domega_c)<eps .OR. abs(domega_r)<eps ) then; rr_fac(iw) = 0d0
     else;                                                rr_fac(iw) = log( domega_r/domega_c )
     endif
     ! rl_fac(his_C(is)) = \int^{his_C}^{his_L} d omega' /(his_C(is) -omega') !  rl_fac(iw) = log( abs((domega_c/domega_l)) )
     if( abs(domega_c)<eps .OR. abs(domega_l)<eps ) then; rl_fac(iw) = 0d0 !
     else ;                                               rl_fac(iw) = log( domega_c/domega_l)
     endif
     if(iw==0) cycle !!  do iw = -nwhis, nwhis !symmetric version. iw=0 is meaningless
     domega_c = zz - his_C(iw)
     if    (iw==nwhis) then;   delta_r = his_R(iw)   - his_C(iw)
     elseif(iw== -1)   then;   delta_r =   0d0       - his_C(iw)
     else;                     delta_r = his_C(iw+1) - his_C(iw)
     endif
     if(iw== -nwhis) then;     delta_l = his_C(iw)  - his_L(iw)
     elseif(iw==  1) then;     delta_l = his_C(iw)  - 0d0
     else;                     delta_l = his_C(iw)  - his_C(iw-1)
     endif
     !          ddr = (his_R(iw)-his_C(iw))/delta_r
     !          ddl = (his_C(iw)-his_L(iw))/delta_l
     rmat(iw)  = rmat(iw  ) + rr_fac(iw)*( 1d0-domega_C/delta_r) !+ ddr
     if(iw/=nwhis .AND. iw/=-1) rmat(iw+1) = rmat(iw+1) + rr_fac(iw)*domega_C/delta_r     !- ddr
     rmat(iw)  = rmat(iw) + rl_fac(iw)*( 1d0+domega_C/delta_l)   !- ddl
     if(iw/=-nwhis .AND. iw/=1) rmat(iw-1) = rmat(iw-1) - rl_fac(iw)*domega_C/delta_l     !+ ddl
  enddo
end subroutine hilbertmat
