!> Calculate W-v zxqi(on the imaginary axis) and zxq(real axis) from sperctum weight rcxq.
module m_dpsion
  public dpsion5
  private
  real(8),allocatable :: his_L(:),his_R(:),his_C(:),rmat(:,:,:),rmatt(:,:,:),rmattx(:,:,:,:),imatt(:,:,:)
  complex(8),allocatable :: imattC(:,:,:)
  real(8),allocatable,save:: gfmat(:,:)
  logical:: eginit=.true.
contains
  subroutine dpsion5(realomega,imagomega,rcxq,nmbas1,nmbas2, zxq,zxqi, chipm,schi,isp,ecut,ecuts) 
    use m_freq,only:  frhis, freqr=>freq_r,freqi=>freq_i, nwhis, npm, nw_i, nw_w=>nw, niwt=>niw
    use m_readgwinput,only: egauss
!    use m_GaussianFilter,only: GaussianFilter
    use m_ftox
    use m_lgunit,only:stdo
    use m_kind,only:kindrcxq
    implicit none
    intent(in)::     realomega,imagomega,     nmbas1,nmbas2,           chipm,schi,isp,ecut,ecuts
    intent(out)::                        rcxq,                zxq,zxqi
    !                                    rcxq is destroyed
    !  works for timereversal=F (npm=2 case).
    !input
    !i   frhis(1:nwhis+1) : specify histgram bins i-th bin is [frhis(i), frhis(i+1)].
    !i   rcxq: the spectrum weight for given bins along the real-axis.
    !i   freqr (0:nw_w) : Calcualte zxq for these real energies.
    !i   freqi (1:niwt) : Calcualte zxqi for these imaginary energies.
    !i   realomega  : A switch to calculate zxq or not.
    !i   imagomega: : A switch to calculate zxqi or not.
    !o   zxq:  W-v along the real axis on freqr(0:nw_w). not accumlating
    !o   zxqi: W-v along the imag axis on freqi(niwt). not accumlating
    !r  We suppose "freqr(i)=moddle of i-th bin; freqr(0)=0." (I think called routine hilbertmat itself is not limited by this condition).
    integer:: igb1,igb2, iw,iwp,ix,ifxx,nmbas1,nmbas2,isp,ispx,it, ii,i,ibas1,ibas2,nmnm
    logical :: evaltest     
    real(8):: px,omp,om,om2,om1, aaa,d_omg, ecut,ecuts,wcut,dee,schi, domega_r,domega_c,domega_l,delta_l,delta_r
    complex(8):: zxq(nmbas1,nmbas2, nw_i:nw_w),zxqi(nmbas1,nmbas2,niwt),img=(0d0,1d0),beta,wfac, zz,rrr(-nwhis:nwhis)
    logical :: realomega, imagomega,chipm,debug=.false.
    integer:: jpm,ipm,verbose,isgi   !     complex(8):: x0mean(nw_i:nw_w,nmbas,nmbas)
    real(8),parameter:: pi  = 4d0*datan(1d0)
    logical::init=.true.
    integer:: imbas1,imbas2
    complex(8):: rcxqin(1:nwhis)
    complex(kindrcxq):: rcxq(nmbas1,nmbas2, nwhis,npm)

    write(stdo,ftox)" -- dpsion5: start... nw_w nwhis=",nw_w,nwhis
    if(chipm.and.npm==2) call rx( 'x0kf_v4h:npm==2 .AND. chipm is not meaningful probably')  ! Note rcxq here is negative 
    call cputid(0)
    GaussianFilter: if(abs(egauss)>1d-15) then
       if(eginit) then
          write(6,'("GaussianFilterX0= ",d13.6)') egauss
          allocate(gfmat(nwhis,nwhis))
          gfmat=gaussianfilterhis(egauss,frhis,nwhis)
          eginit=.false.
       endif
       do ipm=1,npm
          do imbas1=1,nmbas1
             do imbas2=1,nmbas2
                rcxqin = rcxq(imbas1,imbas2,1:nwhis,ipm)
                rcxq(imbas1,imbas2,1:nwhis,ipm) = matmul(gfmat,rcxqin)
             enddo
          enddo
       enddo       !write(6,"(' End of Gaussian Filter egauss=',f9.4)") egauss
    endif GaussianFilter
       
    ispx = merge(isp,3-isp,schi>=0) !  if(schi<0)  ispx = 3-isp  
    if(realomega.and.nwhis <= nw_w) call rxii('dpsion5: nwhis<=nw_w',nwhis,nw_w)
    if(realomega.and.freqr(0)/=0d0) call rx( 'dpsion5: freqr(0)/=0d0') ! I think current version allows any freqr(iw), independent from frhis.

    if(init) then !get Hilbert transformation matrix
      allocate( his_L(-nwhis:nwhis),source=[-frhis(nwhis+1:1+1:-1),0d0,frhis(1  :nwhis)  ])
      allocate( his_R(-nwhis:nwhis),source=[-frhis(nwhis  :1  :-1),0d0,frhis(1+1:nwhis+1)])
      allocate( his_C(-nwhis:nwhis),source=(his_L+his_R)/2d0) !bins are [his_Left,his_Right] !his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
      realomegacase: if(realomega)then;     write(stdo,*) " --- realomega --- "
        if(npm==1) then
          allocate(rmat(0:nw_w,-nwhis:nwhis,npm), source=0d0)
          do it =  0,nw_w
            zz = freqr(it)
            call hilbertmat(zz,  nwhis,his_L,his_C,his_R, rrr)
            rmat(it,:,1) = dreal(rrr)/pi
          enddo
          if(chipm) then
            allocate( rmattx(0:nw_w,nwhis,npm,2) )
            rmattx(:,1:nwhis,1,1) =  rmat(:,1:nwhis,1)
            rmattx(:,1:nwhis,1,2) = -rmat(:,-1:-nwhis:-1,1)
          else  
            allocate( rmatt(0:nw_w,nwhis,npm) )
            rmatt(:,1:nwhis,1) =  rmat(:,1:nwhis,1) - rmat(:,-1:-nwhis:-1,1)
          endif
          deallocate(rmat)
        else  ! npm==2 
          allocate(rmatt(-nw_w:nw_w,nwhis,npm))
          do it  =  -nw_w,nw_w
            zz = merge(-freqr(-it),freqr(it),it<0) 
            call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
            rmatt(it,:,1) =  dreal(rrr  (1:nwhis))/pi
            rmatt(it,:,2) = -dreal(rrr(-1:-nwhis:-1))/pi
          enddo
        endif
      endif realomegacase
      imagomecacase: if(imagomega) then
        if(npm==1) then
          allocate( imatt(niwt, nwhis,npm) )
          do it =  1,niwt
            zz = img*freqi(it)  
            call hilbertmat(zz,nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
            imatt(it,1:nwhis,1) = dreal(rrr(1:nwhis) - rrr(-1:-nwhis:-1))/pi
          enddo
        else ! npm=2 case 
          allocate( imattC(niwt, nwhis,npm) )
          do it =  1,niwt
            zz = img*freqi(it)  
            call hilbertmat(zz,nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
            imattC(it,1:nwhis,1) =   rrr( 1: nwhis   )/pi
            imattC(it,1:nwhis,2) = - rrr(-1:-nwhis:-1)/pi
          enddo
        endif
      endif imagomecacase
      init=.false.
    endif

    do iw= 1, nwhis
      wfac=merge(exp(-(his_C(iw)/ecut)**2 ),1d0, ecut<1d9)     ! rcxq= Average value of Im chi.    Note rcxq is "negative" (
      rcxq(:,:,iw,:)= -wfac/(his_r(iw)-his_l(iw))*rcxq(:,:,iw,:)
    enddo
    if(realomega) then !Hilbert Transformation to get real part
      if(chipm.AND.ispx==1) zxq(:,:,1:nw_w)= zxq(:,:,1:nw_w)+ img*rcxq(:,:,1:nw_w,1) 
      if(.not.chipm)        zxq(:,:,1:nw_w)= img*rcxq(:,:,1:nw_w,1)
      nmnm=2*nmbas1*nmbas2
      if(npm==1.and.chipm) then
        call dgemm('n','t',nmnm,nw_w+1,    nwhis,1d0,dcmplx(rcxq),         nmnm,rmattx(:,:,:,ispx),nw_w+1,1d0,zxq,nmnm)
      elseif(npm==1) then
        call dgemm('n','t',nmnm,nw_w+1,    nwhis,1d0,dcmplx(rcxq),         nmnm,rmatt,             nw_w+1,1d0,zxq,nmnm)
      elseif(npm==2) then
         zxq(:,:,-1:-nw_w:-1)=zxq(:,:,-1:-nw_w:-1) + img*rcxq(:,:,1:nw_w,2)
         !call zaxpy( nmbas1*nmbas2, img, rcxq(1,1,iw,2),1, zxq(:,:,-iw),1)
        call dgemm('n','t',nmnm,npm*nw_w+1,nwhis,1d0,dcmplx(rcxq(:,:,:,1)),nmnm,rmatt(:,:,1),npm*nw_w+1,1d0,zxq,nmnm)
        call dgemm('n','t',nmnm,npm*nw_w+1,nwhis,1d0,dcmplx(rcxq(:,:,:,2)),nmnm,rmatt(:,:,2),npm*nw_w+1,1d0,zxq,nmnm)
      endif
    endif
    if(imagomega) then !Hilbert Transformation to get real part
      nmnm=nmbas1*nmbas2
      if(npm==1) then
        call dgemm('n','t',2*nmnm,niwt,nwhis,1d0,dcmplx(rcxq), 2*nmnm, imatt, niwt, 0d0, zxqi, 2*nmnm )
      elseif(npm==2) then
        call zgemm('n','t', nmnm,niwt,nwhis,1d0, dcmplx(rcxq(:,:,:,1)),nmnm,imattC(1,1,1),niwt, 0d0,zxqi, nmnm )
        call zgemm('n','t', nmnm,niwt,nwhis,1d0, dcmplx(rcxq(:,:,:,2)),nmnm,imattC(1,1,2),niwt, 1d0,zxqi, nmnm )
      endif
    endif
    write(stdo,'("         end dpsion5 ",$)')
    call cputid(0)
  end subroutine dpsion5

!  subroutine GaussianFilter(rcxq,nmbas1,nmbas2, egauss,iprint)
  pure function gaussianfilterhis(egauss, frhis,nwhis) result(gfmat)
    implicit none
    integer,intent(in):: nwhis
    real(8),intent(in):: egauss,frhis(nwhis+1)
    real(8):: gfmat(nwhis,nwhis)
    real(8),allocatable:: frc(:),gfm(:)
    real(8):: ggg
    integer:: i,j
    allocate(frc(nwhis),gfm(nwhis))
    do i=1,nwhis
       frc(i)=(frhis(i)+frhis(i+1))/2d0
    enddo
    do i=1,nwhis
       do j=1,nwhis
          gfm(j)= exp( -(frc(i)-frc(j))**2/(2d0*egauss))
       enddo
       ggg = sum(gfm(:))
       do j=1,nwhis
          gfmat(j,i)= gfm(j)/ggg
       enddo
       
       ! do j=1,nwhis
       !    gfm(j)= frc(j) * exp( -(frc(i)-frc(j))**2/(2d0*egauss))
       ! enddo
       ! ggg = sum(gfm) ! omega*e2(omega) sum rule
       ! do j=1,nwhis
       !    gfmat(j,i) = gfm(j)/frc(j)/ggg
       ! enddo
    enddo
    deallocate(frc,gfm)
  end function gaussianfilterhis
end module m_dpsion

!> Martix for hilbert transformation, rmat.
pure subroutine hilbertmat(zz,nwhis, his_L,his_C,his_R, rmat)
  implicit none
  intent(in):: zz,nwhis,his_L,his_C,his_R
  intent(out):: rmat
  !r  zz is real--->  no img*delta function part
  !r   zz is complex (and Im(zz)>0) : includes all contribution when Im(zz)>eps
  !o  rmat(-nwhis:nwhis) : rmat(0) is not meaningful.
  !i i-th Histgram bin on real axis are given by [his_L, his_R]. center is his_C.
  !r f(zz) = \int_-x(nwhis)^x(nwhis) f(x)/(zz-x)
  !r       = \sum_{i/=0} rmat(i)*f(i)   ,where f(i) is the average value at i-th bin.
  integer :: iw,nwhis
  real(8),parameter:: eps=1d-8, epsz=1d-13
  complex(8),parameter:: img=(0d0,1d0)
  real(8):: his_L(-nwhis:nwhis),his_C(-nwhis:nwhis),his_R(-nwhis:nwhis), delta_r,delta_l,ddr,ddl
  complex(8)::zz,imgepsz, rr_fac(-nwhis:nwhis),rl_fac(-nwhis:nwhis), domega_c,domega_r,domega_l, rmat(-nwhis:nwhis)
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

