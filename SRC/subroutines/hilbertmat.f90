subroutine hilbertmat(zz,nwhis, his_L,his_C,his_R, rmat) ! Martix for hilbert transformation, rmat.
  !r  zz is real--->  no img*delta function part
  !r   zz is complex (and Im(zz)>0) : includes all contribution when Im(zz)>eps
  !o  rmat(-nwhis:nwhis) : rmat(0) is not meaningful.
  !i i-th Histgram bin on real axis are given by [his_L, his_R]. center is his_C.
  !r f(zz) = \int_-x(nwhis)^x(nwhis) f(x)/(zz-x)
  !r       = \sum_{i/=0} rmat(i)*f(i)
  !r     ,where f(i) is the average value at i-th bin.
!!!! 23May2006 I think
!!!! rmat is --------------
!!!! f(zz) = - \int_-x(nwhis)^x(nwhis) f(x)/(zz-x)
!!!!       = - \sum_{i/=0} rmat(i)*f(i)
  ! I forgot minus sign in the previous note.
  !-------------------------------
  implicit none
  integer(4):: iw,nwhis
  complex(8) ::zz,imgepsz
  real(8)    :: his_L(-nwhis:nwhis),his_C(-nwhis:nwhis),his_R(-nwhis:nwhis)
  complex(8) :: rr_fac(-nwhis:nwhis),rl_fac(-nwhis:nwhis),img=(0d0,1d0)
  real(8)::  eps=1d-8, epsz=1d-13,delta_r,delta_l,ddr,ddl
  complex(8):: domega_c,domega_r,domega_l
  complex(8) ::  rmat(-nwhis:nwhis)
  imgepsz =img*epsz
  do iw = -nwhis, nwhis
     if(iw==0) cycle
     domega_r = zz - his_R(iw) + imgepsz
     domega_c = zz - his_C(iw) + imgepsz
     domega_l = zz - his_L(iw) + imgepsz
     if( abs(domega_c)<eps .OR. abs(domega_r)<eps ) then
        rr_fac(iw) = 0d0
     else
        ! rr_fac(his_C(is)) = \int^{his_R}_{his_C} d omega' /(his_C(is) -omega')
        !            rr_fac(iw) = log( abs((domega_r/domega_c)) )
        rr_fac(iw) = log( domega_r/domega_c )
     endif
     if( abs(domega_c)<eps .OR. abs(domega_l)<eps ) then
        rl_fac(iw) = 0d0
     else
        ! rl_fac(his_C(is)) = \int^{his_C}^{his_L} d omega' /(his_C(is) -omega')
        !            rl_fac(iw) = log( abs((domega_c/domega_l)) )
        rl_fac(iw) = log( domega_c/domega_l)
     endif
  enddo
  rmat=0d0
  do iw = -nwhis, nwhis !symmetric version. iw=0 is meaningless
     if(iw==0) cycle
     !          if(debug) print *,' it iw=',it, iw
     domega_c = zz - his_C(iw)
     if(iw==  nwhis) then
        delta_r = his_R(iw)   - his_C(iw)
     elseif(iw== -1) then
        delta_r =   0d0       - his_C(iw)
     else
        delta_r = his_C(iw+1) - his_C(iw)
     endif
     !         if(debug) print *,' it iw RRR1'
     if(iw== -nwhis) then
        delta_l = his_C(iw)  - his_L(iw)
     elseif(iw==  1) then
        delta_l = his_C(iw)  - 0d0
     else
        delta_l = his_C(iw)  - his_C(iw-1)
     endif
     !         if(debug) print *,' it iw RRR2'
     !          ddr = (his_R(iw)-his_C(iw))/delta_r
     !          ddl = (his_C(iw)-his_L(iw))/delta_l
     rmat(iw)  = rmat(iw  ) + rr_fac(iw)*( 1d0-domega_C/delta_r) !+ ddr
     if(iw/=nwhis .AND. iw/=-1) then
        rmat(iw+1) = rmat(iw+1) + rr_fac(iw)*domega_C/delta_r     !- ddr
     endif
     rmat(iw)  = rmat(iw) + rl_fac(iw)*( 1d0+domega_C/delta_l)   !- ddl
     if(iw/=-nwhis .AND. iw/=1) then
        rmat(iw-1) = rmat(iw-1) - rl_fac(iw)*domega_C/delta_l     !+ ddl
     endif
  enddo
end subroutine hilbertmat
