!> Calculate W-v zxqi(on the imaginary axis) and zxq(real axis) from sperctum weight rcxq.
subroutine dpsion5(realomega,imagomega,rcxqin,nmbas1,nmbas2, zxq,zxqi, chipm,schi,isp,ecut,ecuts) 
  use m_freq,only:  frhis, freqr=>freq_r,freqi=>freq_i, nwhis, npm, nw_i, nw_w=>nw, niwt=>niw
  use m_readgwinput,only: egauss
  use m_GaussianFilter,only: GaussianFilter
  use m_lgunit,only:stdo
  implicit none
  intent(in)::     realomega,imagomega,rcxqin,nmbas1,nmbas2,           chipm,schi,isp,ecut,ecuts
  intent(out)::                                              zxq,zxqi
  !     r v4 works for timereversal=F (npm=2 case).
  !     r  See rcxq_zcxq for rcxq, which contains the spectrum weight for given bins along the real-axis.
  !     r ! Note that zxq and zxqi are not accumlating
  !     i frhis(1:nwhis+1) :: specify histgram bins i-th bin is [frhis(i), frhis(i+1)].
  !     i          We suppose "freqr(i)=moddle of i-th bin; freqr(0)=0."
  !     i          (I think called routine hilbertmat itself is not limited by this condition).
  !     i freqr (0:nw_w) : Calcualte zxq for these real energies.
  !     i freqi (1:niwt) : Calcualte zxqi for these imaginary energies.
  !     i   realomega  : A switch to calculate zxq or not.
  !     i   imagomega: : A switch to calculate zxqi or not.
  !     iw rcxq may be altered ---used as work area.
  !     io   zxq :  W-v along the real axis on freqr(0:nw_w)
  !     io   zxqi:  W-v along the imag axis on freqi(niwt)
  !     !
  !     1 Feb2006:  v4 for timereversal=F
  integer:: igb1,igb2, iw,iwp,ix,ifxx,nmbas1,nmbas2,isp,ispx,it, ii,i,ibas1,ibas2,nmnm
  logical :: evaltest     
  real(8):: px,omp,om,om2,om1, aaa,d_omg, ecut,ecuts,wcut,dee,schi, domega_r,domega_c,domega_l,delta_l,delta_r
  complex(8):: rcxq(nmbas1,nmbas2, nwhis,npm) ,rcxqin(nmbas1,nmbas2, nwhis,npm) 
  complex(8):: zxq(nmbas1,nmbas2, nw_i:nw_w),zxqi(nmbas1,nmbas2,niwt),img=(0d0,1d0),beta,wfac, zz
  logical :: realomega, imagomega,chipm,debug=.false.
  real(8),allocatable :: his_L(:),his_R(:),his_C(:),rmat(:,:,:),rmati(:,:,:),rmatt(:,:,:),imatt(:,:,:)
  complex(8),allocatable :: rmatiC(:,:,:),imattC(:,:,:),zxqn(:,:),zxqn1(:,:,:),rx0mean1(:,:,:),rx0mean(:), rrr(:)
  integer:: jpm,ipm,verbose,isgi   !     complex(8):: x0mean(nw_i:nw_w,nmbas,nmbas)
  real(8),allocatable:: ebb(:)
  real(8),parameter:: pi  = 4d0*datan(1d0)
  if(chipm.and.npm==2) call rx( 'x0kf_v4h:npm==2 .AND. chipm is not meaningful probably')  ! Note rcxq here is negative 
  write(stdo,'(" -- dpsion5: start...  nw_w nwhis=",2i5)') nw_w,nwhis
  call cputid(0)
  if(abs(egauss)>1d-15) call GaussianFilter(rcxq,nmbas1,nmbas2,egauss,iprint=.true.) !Smearging Imag(X0). Use egauss = 0.05 a.u.\sim 1eV for example.
  ispx = merge(isp,3-isp,schi>=0) !  if(schi<0)  ispx = 3-isp  
  if(realomega.and.nwhis <= nw_w) call rxii( ' dpsion5: nwhis<=nw_w',nwhis,nw_w)
  if(realomega.and.freqr(0)/=0d0) call rx( ' dpsion5: freqr(0)/=0d0') ! I think current version allows any freqr(iw), independent from frhis.
  if(debug) write(stdo,*)' dpsion5: RRR 2222222222 sumchk 111 rcxq=', sum(abs(rcxq))
  allocate( his_L(-nwhis:nwhis),source=[-frhis(nwhis+1:1+1:-1),0d0,frhis(1  :nwhis)  ])
  allocate( his_R(-nwhis:nwhis),source=[-frhis(nwhis  :1  :-1),0d0,frhis(1+1:nwhis+1)])
  allocate( his_C(-nwhis:nwhis),source=(his_L+his_R)/2d0) ! bins are [his_Left,his_Right] !his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
  rcxq=rcxqin
  do iw= 1, nwhis
     wfac=merge(exp(-(his_C(iw)/ecut)**2 ),1d0, ecut<1d9)  ! rcxq is used as work---> rcxq= Average value of Im chi.    Note rcxq is "negative" (
     rcxq(:,:,iw,:)= -wfac/(his_r(iw)-his_l(iw))*rcxq(:,:,iw,:)
  enddo
  if(debug) write(stdo,*)'sumchk 122 rcxq=', sum(abs(rcxq))
  ! debugevaltest: if(evaltest() .AND. nmbas1==nmbas2) then; write(stdo,"('hhh --- EigenValues for rcxq --------')")
  !    allocate(ebb(nmbas1))
  !    do jpm= 1,npm; do iw = 1, nwhis;           call diagcvh2(rcxq(:,:,iw,jpm),nmbas1,ebb)
  !       do ii=1,nmbas1
  !          write(stdo,"('hhh1: xxxxxxxxxxxxxxxxx',2i4)") jpm,iw
  !          if(abs(ebb(ii))>1d-8 .AND. ebb(ii)>0) write(stdo,"('hhh1: jpm iw eb=',2i4,d13.5)") jpm,iw,ebb(ii)
  !       enddo
  !    enddo;  enddo
  !    deallocate(ebb)
  ! endif debugevaltest
  realomegacase: if(realomega)then;     write(stdo,*) " --- realomega --- "
     if(npm==1) then
        allocate(rmat(0:nw_w,-nwhis:nwhis,npm), source=0d0)
        allocate(rrr(-nwhis:nwhis),source=(0d0,0d0))
        do it =  0,nw_w
           zz = freqr(it)
           call hilbertmat(zz,  nwhis,his_L,his_C,his_R, rrr)
           rmat(it,:,1) = dreal(rrr)
        enddo
        allocate( rmatt(0:nw_w,nwhis,npm) )
        if( chipm .AND. ispx==1 ) rmatt(:,1:nwhis,1) =  rmat(:,1:nwhis,1)
        if( chipm .AND. ispx==2 ) rmatt(:,1:nwhis,1) = -rmat(:,-1:-nwhis:-1,1)
        if(.not.chipm)            rmatt(:,1:nwhis,1) =  rmat(:,1:nwhis,1) - rmat(:,-1:-nwhis:-1,1)
        deallocate(rmat,rrr)
     else  ! npm==2 
        allocate(rmatt(-nw_w:nw_w,nwhis,npm))
        allocate(rrr(-nwhis:nwhis),source=(0d0,0d0))
        do it  =  -nw_w,nw_w
           zz = merge(-freqr(-it),freqr(it),it<0) 
           call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
           rmatt(it,:,1) =  dreal(rrr  (1:nwhis))
           rmatt(it,:,2) = -dreal(rrr(-1:-nwhis:-1))
        enddo
        deallocate(rrr)
     endif
     rmatt = rmatt/pi ; if(debug) write(stdo,*)'dpsion5: RRR 6666'  ! WARN! I think npm==2.and.chipm does not make sense. apr2012.
     if(chipm.AND.ispx==1) zxq(:,:,1:nw_w)= zxq(:,:,1:nw_w)+ img*rcxq(:,:,1:nw_w,1) 
     if(.not.chipm)        zxq(:,:,1:nw_w)= img*rcxq(:,:,1:nw_w,1)
     nmnm=2*nmbas1*nmbas2
     if(npm==1) call dgemm('n','t',nmnm,nw_w+1,    nwhis,1d0,rcxq,         nmnm,rmatt,           nw_w+1,1d0,zxq,nmnm)
     if(npm==2) then
        zxq(:,:,-1:-nw_w:-1)=zxq(:,:,-1:-nw_w:-1) + img*rcxq(:,:,1:nw_w,2) !call zaxpy( nmbas1*nmbas2, img, rcxq(1,1,iw,2),1, zxq(:,:,-iw),1)
        call dgemm('n','t',nmnm,npm*nw_w+1,nwhis,1d0,rcxq(1,1,1,1),nmnm,rmatt(:,:,1),npm*nw_w+1,1d0,zxq,nmnm)
        call dgemm('n','t',nmnm,npm*nw_w+1,nwhis,1d0,rcxq(1,1,1,2),nmnm,rmatt(:,:,2),npm*nw_w+1,1d0,zxq,nmnm)
     endif   
     if(.not.(npm==1.or.npm==2)) call rx( 'dpsion5: npm=1 or 2')
     deallocate(rmatt)
  endif realomegacase
  imagomecacase: if(imagomega) then
     nmnm=nmbas1*nmbas2
     allocate(rrr(-nwhis:nwhis))
     if(npm==1) then
        allocate( rmati (niwt,-nwhis:nwhis,npm),source=0d0)
        do it =  1,niwt
           zz = img*freqi(it)  
           call hilbertmat(zz,nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
           rmati (it,:,1) = dreal(rrr)
        enddo
        allocate( imatt(niwt, nwhis,npm) )
        imatt(:,1:nwhis,1) = (rmati(:,1:nwhis,1) - rmati(:,-1:-nwhis:-1,1))/pi
        deallocate(rmati,rrr)
        call dgemm('n','t',  2*nmnm, niwt, nwhis, 1d0, rcxq, 2*nmnm, imatt, niwt, 0d0, zxqi, 2*nmnm )
        deallocate(imatt)
     else ! npm=2 case 
        allocate( rmatiC(niwt,-nwhis:nwhis,npm),source=(0d0,0d0))
        do it =  1,niwt
           zz = img*freqi(it)  
           call hilbertmat(zz,nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
           rmatiC(it,:,1) = rrr
        enddo
        allocate( imattC(niwt, nwhis,npm) )
        imattC(:,1:nwhis,1) =   rmatiC(:, 1: nwhis,   1)/pi
        imattC(:,1:nwhis,2) = - rmatiC(:,-1:-nwhis:-1,1)/pi
        call zgemm('n','t', nmnm,niwt,nwhis,1d0,rcxq(1,1,1,1),nmnm,imattC(1,1,1),niwt, 0d0,zxqi, nmnm )
        call zgemm('n','t', nmnm,niwt,nwhis,1d0,rcxq(1,1,1,2),nmnm,imattC(1,1,2),niwt, 1d0,zxqi, nmnm )
        deallocate(rmatiC,rrr,imattC)
     endif
  endif imagomecacase
  deallocate(his_L,his_C,his_R)
  write(stdo,'("         end dpsion5 ",$)')
  call cputid(0)
end subroutine dpsion5
