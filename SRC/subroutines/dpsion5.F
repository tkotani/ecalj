      subroutine dpsion5(       !frhis,nwhis, freqr,nw_w, freqi,niwt, 
     i     realomega,   imagomega, rcxq, nmbas1,nmbas2, !npm, nw_i,
     o     zxq,zxqi,
     i     chipm,schi,isp, ecut,ecuts)
      use m_freq,only:  frhis, freqr=>freq_r,freqi=>freq_i, nwhis, npm, nw_i, nw_w=>nw, niwt=>niw
      use m_readgwinput,only: egauss
      use m_GaussianFilter,only: GaussianFilter
      implicit none
      intent(in)::
     i     realomega,   imagomega, rcxq, nmbas1,nmbas2, 
     i     chipm,schi,isp, ecut,ecuts
      intent(out)::
     o     zxq,zxqi
C     - Calculate W-v zxqi(on the imaginary axis) and zxq(real axis) from sperctum weight rcxq.
C     r v4 works for timereversal=F (npm=2 case).
C     r  See rcxq_zcxq for rcxq, which contains the spectrum weight for given bins along the real-axis.
C     r ! Note that zxq and zxqi are not accumlating
C     i frhis(1:nwhis+1) :: specify histgram bins i-th bin is [frhis(i), frhis(i+1)].
C     i          We suppose "freqr(i)=moddle of i-th bin; freqr(0)=0."
C     i          (I think called routine hilbertmat itself is not limited by this condition).
C     i freqr (0:nw_w) : Calcualte zxq for these real energies.
C     i freqi (1:niwt) : Calcualte zxqi for these imaginary energies.
C     i   realomega  : A switch to calculate zxq or not.
C     i   imagomega: : A switch to calculate zxqi or not.
C     iw rcxq may be altered ---used as work area.
C     io   zxq :  W-v along the real axis on freqr(0:nw_w)
C     io   zxqi:  W-v along the imag axis on freqi(niwt)
C     !
C     1 Feb2006:  v4 for timereversal=F
C     ! July2005: v3Add spin chipm mode
C     ! July2005: This version alter rcxq----it is used as work area.
C     ! sergey faleev Apr 2002 ; Rebuiled by takao
C------------------------------------------------------------------
      integer(4):: igb1,igb2, iw,iwp,ix,ifxx,nmbas1,nmbas2 ! nw_w,niwt,,nwhis,npm,nw_i,
      real(8) :: pi,px,omp,om,om2,om1, 
     &     aaa,d_omg            !frhis(nwhis+1), freqr(0:nw_w),  freqi(niwt),
      logical :: realomega, imagomega
      complex(8):: rcxq(nmbas1,nmbas2, nwhis,npm) !sf 13June
c     logical   :: iepsmode
      logical :: chipm
      integer(4)::isp,ispx      !, nmbas
c     complex(8):: rcxqmean(nwhis,npm,nmbas,nmbas)  !takao sep2006 add nmbas
C...  ecut mode
      real(8):: ecut,ecuts,wcut,wcutef,dee,schi
      logical ::debug=.false.
      real(8),allocatable :: his_L(:),his_R(:),his_C(:)
      integer(4) it
      real(8):: domega_r,domega_c,domega_l,delta_l,delta_r
      real(8),allocatable ::rmat(:,:,:),rmati(:,:,:),rmatt(:,:,:),imatt(:,:,:)
      complex(8),allocatable :: rmatiC(:,:,:),imattC(:,:,:)
      complex(8) ::beta,wfac
      complex(8):: zz
      complex(8),allocatable :: zxqn(:,:),zxqn1(:,:,:),rx0mean1(:,:,:),rx0mean(:)
      complex(8),allocatable:: rrr(:)
      integer(4):: jpm,ipm,verbose,isgi
c     complex(8):: x0mean(nw_i:nw_w,nmbas,nmbas)
      complex(8):: 
     o     zxq (nmbas1,nmbas2, nw_i: nw_w), !iw=0 means omg=0,
!iw=1:nw_w corresponds to iw's bit of the frequensy histogram 
     o     zxqi(nmbas1,nmbas2,niwt),img !npm), img  !zxqi(...,npm) may2006
      real(8),allocatable:: ebb(:)
      integer(4):: ii,i,ibas1,ibas2
      logical :: evaltest       !,testtr
c     ------------------------------------------------
      write(6,'(" -- dpsion5: start...   ",$)')
      write(6,"('  nw_w nwhis=',2i5)") nw_w,nwhis
      
!     ! Gaussian filtering of rcxq. Smearging Imag(X0). We may use egauss = 0.05 a.u.\sim 1eV for example.
!     ! Jun2020
      if(abs(egauss)>1d-15) then
         call GaussianFilter(rcxq,nmbas1,nmbas2,egauss,iprint=.true.)
      endif   

      if(debug) then
         write(6,*)' nmbas1 nmbas2 nwhis npm =',  nmbas1,nmbas2,nwhis,npm
         write(6,*)' sumchk rcxq=', sum(abs(rcxq))
      endif
      pi  = 4d0*datan(1d0)
      img = (0d0,1d0)
      call cputid(0)
      ispx = isp
      if(schi<0) then
         ispx = 3-isp           !flip
      endif

!     ! Check freqr
      if(realomega) then
         if( nwhis <= nw_w ) then
            write(6,*)nwhis,nw_w
            call rx( ' dpsion5: nwhis<=nw_w')
         endif
         if( freqr(0)/=0d0 ) call rx( ' dpsion5: freqr(0)/=0d0')
!     ! I think current version allows any freqr(iw), independent from frhis.
c$$$  aaa = 0d0
c$$$  if(nw_w>0) then
c$$$  do iw = 1,nw_w
c$$$  aaa = aaa + abs( freqr(iw) - (frhis(iw)+frhis(iw+1))/2d0 )
c$$$  if(debug) write(6,"(' iw freqr frhis_m=',i5,2f13.6)" )
c$$$  &        iw,freqr(iw),  (frhis(iw)+frhis(iw+1))/2d0
c$$$  enddo
c$$$  if(aaa>1d-10)call rx( 'dpsion5:freqr/=frhis_m is not implimented yet')
c$$$  endif
      endif                     !realomega

C--------------------------------------------------------------
!     ! Each histogram bins are  [his_Left, his_Right], and  his_Center is middle.
!     ! his_C(0) is at zero. his_R(0) and his_L(0) are not defined.
      if(debug) write(6,*)' dpsion5: RRR 2222222222 '
      allocate(his_L(-nwhis:nwhis),his_R(-nwhis:nwhis),his_C(-nwhis:nwhis))
      his_L(1:nwhis) = frhis(  1:  nwhis)
      his_R(1:nwhis) = frhis(1+1:1+nwhis)
      his_C(1:nwhis) = (his_L(1:nwhis) + his_R(1:nwhis) )/2d0
      do iw= 1,nwhis
         his_L(-iw) = -his_R(iw)
         his_R(-iw) = -his_L(iw)
         his_C(-iw) = -his_C(iw)
      enddo
      his_C(0) = 0d0; his_R(0)=-999; his_L(0)=-999
C     
      if(debug) write(6,*)'sumchk 111 rcxq=', sum(abs(rcxq))

      do iw= 1, nwhis
         if(ecut<1d9) then
            wfac= wcutef(his_C(iw), ecut,ecuts)
         else
            wfac= 1d0
         endif
!     rcxq is used as work---> rcxq= Average value of Im chi.
!     Note rcxq is "negative" (
         do jpm=1,npm
            call dscal(2*nmbas1*nmbas2, -wfac/(his_r(iw)-his_l(iw)),rcxq(1,1,iw,jpm),1)
         enddo
c     if(debug) write(6,*) 'dpsion5: RRR 7777 iw wfac=',iw,wfac,ecut,ecuts
      enddo
      if(debug) write(6,*)'sumchk 122 rcxq=', sum(abs(rcxq))
      if(debug) write(6,*)'sumchk 222 rcxq=', sum(abs(rcxq))
      if(evaltest().and.nmbas1==nmbas2) then
         write(6,"('hhh --- EigenValues for rcxq --------')")
         allocate(ebb(nmbas1))
         do jpm= 1,npm
            do iw = 1, nwhis
               call diagcvh2(rcxq(:,:,iw,jpm),nmbas1,ebb)
               do ii=1,nmbas1
                  write(6,"('hhh1: xxxxxxxxxxxxxxxxx',2i4)") jpm,iw
                  if(abs(ebb(ii))>1d-8.and.ebb(ii)>0) 
     &                 write(6,"('hhh1: jpm iw eb=',2i4,d13.5)") jpm,iw,ebb(ii)
               enddo
            enddo
         enddo
         deallocate(ebb)
      endif

C---  realomega case
      if(realomega)then
         write(6,*) " --- realomega --- "
         if(npm==1) then
            allocate( rmat(0:nw_w,-nwhis:nwhis,npm), rrr(-nwhis:nwhis))
            rmat  = 0d0
            do it =  0,nw_w
               zz = freqr(it)   !his_C(it)
               call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
               rmat(it,:,1) = dreal(rrr)
            enddo ;   if(debug) write(6,*) 'dpsion5: RRR 55555555555'
            allocate( rmatt(0:nw_w,nwhis,npm) )
            if(     chipm.and.ispx==1 ) then
               rmatt(:,:,1) = rmat(:,1:nwhis,1)
            elseif( chipm.and.ispx==2 ) then
               do iw= 1,nwhis
                  rmatt(:,iw,1) = -rmat(:,-iw,1)
               enddo
            else
               do iw= 1,nwhis
                  rmatt(:,iw,1) = rmat(:,iw,1) - rmat(:,-iw,1)
               enddo
            endif
            deallocate(rmat,rrr)
         else                   ! npm==2 case -------------------------------------------------
            allocate( rmatt(-nw_w:nw_w,nwhis,npm), rrr(-nwhis:nwhis))
            rmatt = 0d0
            do it  =  -nw_w,nw_w
               if(it<0) then
                  zz = -freqr(-it) !his_C(it)
               else
                  zz = freqr(it) !his_C(it)
               endif
               call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr)
               rmatt(it,:,1) =  dreal(rrr  (1:nwhis))
               rmatt(it,:,2) = -dreal(rrr(-1:-nwhis:-1))
            enddo ;   if(debug) write(6,*) 'dpsion5: RRR2 55555555555'
            deallocate(rrr)
         endif
         rmatt = rmatt/pi ; if(debug) write(6,*)'dpsion5: RRR 6666'
!     ! takao remove if(nolfc) block here.
c     write(6,*) " --- realomega dgemm--- "
!     ! WARN! I think npm==2.and.chipm does not make sense. apr2012.
         if(npm==2.and.chipm) call rx( 'x0kf_v4h:npm==2.and.chipm is not meaningful probably')
!     ! Note rcxq is negative now (converted at the top of this routine !!!
         if(     chipm .and. ispx==2 ) then
!nothing here
!Since the range of zxq is nw_i=0, we have no area to store negative energy part of chipm.
         elseif( chipm               ) then
            call zaxpy( nmbas1*nmbas2*nw_w, img, rcxq, 1, zxq(1,1,1), 1)
         else
            zxq = 0d0           ! not accumlating case.
            call zaxpy( nmbas1*nmbas2*nw_w, img, rcxq(1,1,1,1), 1, zxq(1,1,1), 1)
         endif

         if(npm==2) then
            do iw=1,nw_w
               call zaxpy( nmbas1*nmbas2, img, rcxq(1,1,iw,2),1, zxq(:,:,-iw),1)
            enddo
         endif

         if(npm==1) then
            call dgemm('n','t',  2*nmbas1*nmbas2, nw_w+1, nwhis, 1d0, 
     &           rcxq, 2*nmbas1*nmbas2,  rmatt, nw_w+1,
     &           1d0, zxq, 2*nmbas1*nmbas2 )
         elseif(npm==2) then
            call dgemm('n','t',  2*nmbas1*nmbas2,   npm*nw_w+1, nwhis, 1d0, 
     &           rcxq(1,1,1,1), 2*nmbas1*nmbas2, rmatt(:,:,1), npm*nw_w+1,
     &           1d0, zxq, 2*nmbas1*nmbas2 )
            call dgemm('n','t',  2*nmbas1*nmbas2,   npm*nw_w+1, nwhis, 1d0, 
     &           rcxq(1,1,1,2), 2*nmbas1*nmbas2, rmatt(:,:,2), npm*nw_w+1,
     &           1d0, zxq, 2*nmbas1*nmbas2 )
         else
C     stop2rx 2013.08.09 kino            stop 'dpsion5: npm=1 or 2'
            call rx( 'dpsion5: npm=1 or 2')
         endif
         deallocate(rmatt)
      endif

!     ! === imagomega case      imatt(niwt -->niwt,npm may2005 ===
      if(imagomega) then
         allocate( rrr(-nwhis:nwhis))
         if(npm==1) then
            allocate( rmati (niwt,-nwhis:nwhis,npm))
            rmati= 0d0
         else
            allocate( rmatiC(niwt,-nwhis:nwhis,npm))
            rmatiC = 0d0
         endif ;   if(debug) write(6,*) 'dpsion5: III 111111155555555555'
         do it =  1,niwt
            zz = img*freqi(it)  !his_C(it)
            call hilbertmat(zz, nwhis,his_L,his_C,his_R, rrr) !Im(zz)>0
            if(npm==1) then
               rmati (it,:,1) = dreal(rrr)
            else
               rmatiC(it,:,1) = rrr
            endif
         enddo ;   if(debug) write(6,*) 'dpsion5: III 55555555555'
!     ! ==== npm=1 case ====
         if(npm==1) then
            allocate( imatt(niwt, nwhis,npm) )
            do iw= 1,nwhis
               imatt(:,iw,1) = rmati(:,iw,1) - rmati(:,-iw,1)
            enddo
            deallocate(rmati,rrr)
            imatt = imatt/pi; if(debug) write(6,*) 'dpsion5: III  '
            call dgemm('n','t',  2*nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &           rcxq, 2*nmbas1*nmbas2, imatt, niwt,
     &           0d0, zxqi, 2*nmbas1*nmbas2 )
            deallocate(imatt)
!     ! ==== npm=2 case ====
         else
            allocate( imattC(niwt, nwhis,npm) )
            do iw= 1,nwhis
               imattC(:,iw,1) =   rmatiC(:, iw,1)
               imattC(:,iw,2) = - rmatiC(:,-iw,1)
            enddo
            deallocate(rmatiC,rrr)
            imattC = imattC/pi; if(debug) write(6,*) 'dpsion5: IIIc '
            call zgemm('n','t',  nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &           rcxq(1,1,1,1), nmbas1*nmbas2, imattC(1,1,1), niwt,
     &           0d0, zxqi,    nmbas1*nmbas2 )
            call zgemm('n','t',  nmbas1*nmbas2, niwt, nwhis, 1d0, 
     &           rcxq(1,1,1,2), nmbas1*nmbas2, imattC(1,1,2), niwt,
     &           1d0, zxqi,    nmbas1*nmbas2 )
            deallocate(imattC)
         endif
      endif
      deallocate(his_L,his_C,his_R)
      write(6,'("         end dpsion5 ",$)')
      call cputid(0)
      end
      
      real(8) function wcutef(e,ecut,ecuts)
      real(8):: e,ecut,ecuts
c     wcutef = 1d0/( exp((e-ecut)/ecuts)+ 1d0)
      wcutef = exp( -(e/ecut)**2 ) ! ecuts is not used in this case
      end
