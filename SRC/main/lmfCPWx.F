      program lmfCPW
      implicit none
      integer:: i,j,ikp,ib1,ib2,it,nmx,nev,jsp,nev0
      complex(8)::img=(0d0,1d0),phase,phase0,pp1,pp2,pp3,phase1,phase2,phase3,ooi,ooj
      complex(8),allocatable:: hamm(:,:),ovlm(:,:),t_zv(:,:),fac(:,:),ovlmi(:,:),bb(:),bbx(:),
     & eemat(:,:),bqr(:,:),ovlm0(:,:)
      real(8),allocatable:: evl(:),GG(:,:),ggn(:,:),bqx(:),xqpp(:),qpp(:,:),qsyms(:,:),qsyme(:,:),
     & eee(:,:) ,eeex(:),percent(:,:)
      real(8)::qp(3),pi=4d0*atan(1d0),norm,evalx,alpha,gmax,plat(3,3),qlat(3,3),alat,ggg(3),fac1,fac2,deltar(3),gmx
      logical:: lprint=.true.,savez=.false.,getz=.false. !dummy
      integer:: ifig=-999       !dummy
      integer:: ndatx,ifsy1,ifsy2,ifile_handle,ifsy,ig
      real(8)::Grange,Gscale,gnmax,ggabs,gtestmx,wsum,xx(3),alfac,dx,dqp(3),xqp,efac,dummy,wsr,vol,bohr,scaling,ggamma
      integer:: iev,jev,ixx,nevout,idat,init,iend,isym,nsym,ifiei,ifieg,ifidi,nx
      integer,allocatable:: ndiv(:),ikpe(:)
      integer,allocatable:: igv2x0(:,:)
      integer:: cases
      integer:: napw0,nnn,iapw,nnnn,napw,ix,iy,iz,ifp,m,k,ifx,ikpix(2,100),ikk,ifp2,ifp3,nkpe,ikx,ifp4
      real(8):: eigov(3),diffe,ediff(3)
      character(100)::fff1,fff2,fff3,fname,ftag,fff22
      character(300)::aaa
!!! Settings -------
      bohr= 0.529177
      alat= 5.43095/bohr !for Si
      plat(:,1)= [0.000000,   0.500000,   0.500000]
      plat(:,2)= [0.500000,   0.000000,   0.500000]
      plat(:,3)= [0.500000,   0.500000,   0.000000]
!!!
      ifx=ifile_handle()
      open(ifx,file='ikpix.input')
      do i=1,100
         read(ifx,*,end=111)ikpix(1,i),ikpix(2,i)
         ikk=i
      enddo
 111  continue
      close(ifx)

c     do i=1,ikk
c         write(6,*)ikpix(1,i),ikpix(2,i)
c      enddo
c      stop 
!!!   Get qlat and WignerSeitzRadius
      call dinv33(plat,1,qlat,vol) !vol primitive cell volume in real space.
      wsr = (vol*3d0/(4d0*pi))**(1d0/3d0) !4pi/3 wsr**3 = vol

c      gmax  = 4.1 !4 ==>59 basis 4.1==>65, 6==>229  8==> in unit of 2pi/a
c      gnmax = gmax**2 
c      fac1= 2d0 !Galpha  =      fac1*wsr    ! Galpha \propto number of G domain
c      fac2= 2d0 !Gbeta   = 1d0/(fac2*wsr)   ! 1d0 !1d0  ! Scaling factor to set ggn best result
c      cases = 1  !case for scaling
      nnn=12                   ! for integer lattice covering large enough G vectors.

!!!   Readin factors

      write(6,*)' Grange Gscale Gmax cases=?'
      read(5,*) Grange,Gscale,gmx, cases
c      Gscale   =    ! fac2*1d0/wsr  ! 1d0 !1d0  ! Scaling factor to set ggn best result
c      Grange   =    ! G range wsr
      gnmax = gmx**2 !4.1==>65, 4.5**2, 6**2==>229  in unit of 2pi/a
!!
      
      
!!!   igv2x0 is the integer lattice
      nnnn=(2*nnn+1)**3
      napw=nnnn
      allocate(igv2x0(3,napw))
      iapw=0
      do ix=-nnn,nnn
      do iy=-nnn,nnn
      do iz=-nnn,nnn
         ggg= matmul(qlat,[ix,iy,iz])
         if( sum(ggg**2) < ( gmx + 5*Grange )**2) then
            iapw=iapw+1
            igv2x0(:,iapw)=[ix,iy,iz]
         endif   
c     print *,' iPW=',iapw,igv2x0(:,iapw)
      enddo
      enddo
      enddo
      napw0=iapw
      write(6,"('napw0=',i5)") napw0

     
!!!   Info shown
c     write(6,"(' Gnmax Galpha Gbeta=noscaling',3f13.5)") gnmax,Galpha!,Gbeta
      write(6,"(' Gmax Grange Gscaling=',3f13.5)") gmx,Grange,Gscale
      allocate(evl(napw0),GG(3,napw0), ggn(3,napw0) )
      iev=0
      nx=0
      efac=2 !up to what energy. Empty dispersion.
      do i=1,napw0
         GG(:,i)=  matmul(qlat,igv2x0(:,i)) !G vectors are generated
         if(sum(GG(:,i)**2)< gnmax*efac) then
            nx=nx+1 !nx count is just for empty sphere dispersion
         endif
c         write(6,'(i5,2x,3i3,2x,3f9.3)')i,igv2x0(:,i), matmul(qlat,igv2x0(:,i))

ccccccccccccc
c used with 'add extra base, below
c     if(sum(GG(:,i)**2)< 4.1 .and. sum(GG(:,i)**2)> 1d-6) cycle
ccccccccccccc            
         
         if(sum(GG(:,i)**2)< gnmax) then !Gn generatos
            iev = iev+1
            ggabs = sum(GG(:,i)**2)**.5
            if(cases==0) then
               scaling = 1d0/Gscale
            elseif(cases==1) then
               scaling = (ggabs)/(ggabs+Gscale) 
            elseif(cases==100) then
               scaling = ggabs**0.5/(ggabs**0.5 + Gscale) 
c               scaling = 1d0/( exp(-(ggabs-3d0)) + 1d0 ) !band reasonable but .5d-11
c               scaling = 1d0/( exp(-(ggabs-3d0)*0.5d0) + 1d0 ) ! .1d-11
c               scaling = 1d0/( exp(-(ggabs-3d0)*0.2d0) + 1d0 ) ! .1d-11 !all are around 0.5

               
c               if(ggabs<3d0) then !not good enough band above 20eV
c                  scaling = 0.5d0
c               else
c                  scaling = 1d0
c               endif   

c               scaling = 0.1*exp(- ggabs/Galpha) + (1d0-Galpha/(ggabs+1d-8)) * (1- exp(-ggabs/Galpha)) 
c               scaling = (1d0 - Galpha/(ggabs+1d-8)) * (1 - exp(-0.5*ggabs/Galpha)) 

c               scaling = 1d0 !(ggabs)/(ggabs+0.1) 

c     scaling = (ggabs**2)/(ggabs**2+Gbeta)
            else
               stop 'error not fit to allowed cases'
            endif   
            ggn(:,iev)= GG(:,i)*scaling
            write(6,'("scaling Gn G",i5,i7,x,f8.3,x,3f9.3,x,3f9.3)')iev,i, scaling, ggn(:,iev),GG(:,i)
         endif   
      enddo
      nev=iev
      nev0=nev
      write(6,"('napw0 nev  gnmax Grange Gscale cases=',2i5,3f9.4,l)") napw0,nev, gnmax, Grange, Gscale, cases

c$$$c$$$c$$$c$$$ccccccccccccccccccccccccccccccccccccccccc
c$$$!!   add extra base ---> not improve so much. To get good band, we have eigovratio \sim 1d-6.
c$$$      write(6,*)
c$$$      do i=1,napw0
c$$$         GG(:,i)=  matmul(qlat,igv2x0(:,i)) !G vectors are generated
c$$$         if(sum(GG(:,i)**2)< 1d-6) cycle
c$$$
c$$$         if(sum(GG(:,i)**2)< 4.2) then !Gn generatos
c$$$            iev = iev+1
c$$$            scaling = 0.3
c$$$c            ggabs = sum(GG(:,i)**2)**.5
c$$$c            scaling = ggabs/(ggabs+Gscale) 
c$$$            ggn(:,iev)= GG(:,i)*scaling
c$$$            write(6,'("scaling Gn G",i5,i7,x,f8.3,x,3f9.3,x,3f9.3)')iev,i, scaling, ggn(:,iev),GG(:,i)
c$$$         endif   
c$$$      enddo
c$$$      nev=iev
c$$$      write(6,"('add extra bases: nev=',i5)") nev
c$$$c$$$c$$$c$$$cccccccccccccccccccccccccccccccccccccccccccc      
      
      
!!! --- plat, qlat, vol, wsr(one atom per cell) are to a file (just info).
c      print *,'vol asize=',vol,wsr
      ifp=ifile_handle()
      write(fff1,"(f9.3)")Grange
      write(fff2,"(f9.3)")Gscale
      write(fff22,"(f9.3)")gmx
      write(fff3,"(i3)")cases
      ftag=trim(adjustl(fff1))//'_'//trim(adjustl(fff2))//'_'//trim(adjustl(fff22))//'_'//trim(adjustl(fff3))
      fname='CPWdata.'//trim(ftag)
      open(ifp,file=trim(fname))
      write(ifp,"(3f11.6,5x,3f11.6)") ((plat(m,k),m=1,3),(qlat(m,k),m=1,3),k=1,3)
      write(ifp,"('             PLAT              and         QLAT     ')")
      write(ifp,"(3f12.5,'     ! alat, CellVol, WSradius ')") alat, vol,wsr
      write(ifp,"(f9.4,i9,i9,' ! gmx nev napw0 ')")     gmx,nev,napw0
      write(ifp,"(2f9.4,i5, '  ! Grange Gscale cases')") Grange,Gscale,cases
      write(ifp,"(3f9.4, ' ! Scaling at Gmax, Grange')") gmx/(gmx+Gscale), Grange
c      close(ifp)
      
      
!!!! --- q-point generator: xqpp and qpp for i=1,ndatx
      nsym=2
      allocate(ndiv(nsym))
      ndiv(1)=40
      ndiv(2)=40
      allocate(xqpp(sum(ndiv)+nsym),qpp(3,sum(ndiv)+nsym),qsyms(3,nsym),qsyme(3,nsym),ikpe(0:nsym))
      qsyms(:,1) = [1d0,0d0,0d0]
      qsyme(:,1) = [0d0,0d0,0d0]
      qsyms(:,2) = [0d0,0d0,0d0]
      qsyme(:,2) = [0.5, -0.5, -0.5]
      idat=0
      nkpe=0
      ikpe(0)=1
      do isym= 1,nsym
         dx  = sum( (qsyme(:,isym)-qsyms(:,isym))**2 )**0.5/ndiv(isym)
         dqp = (qsyme(:,isym)-qsyms(:,isym))/ndiv(isym) 
         if(isym==1) then
            qp = qsyms(:,1) - dqp
            xqp= - dx
         endif
         print *,'dx,dp=',dx,dqp
         do ikp= 0,ndiv(isym)
            if(isym>1.and.ikp==0) cycle
            qp  = qp +  dqp
            xqp = xqp + dx
            idat =idat+1
            xqpp(idat)  = xqp
            qpp(:,idat) = qp
            if(ikp == ndiv(isym)) then
               nkpe=nkpe+1
               ikpe(nkpe)=idat
            endif
            print *,' xxxxxx',idat,xqp, qpp(:,idat)
         enddo
      enddo
      ndatx=idat
      do i=1,ndatx
        write(6,"( '  xqp qp=',f10.3,2x,3f10.4)") xqpp(i), qpp(:,i)
      enddo
c      stop 'xxxxxxxx'
!     ! --- open files
      ifieg=ifile_handle()
      open(ifieg,file='EmptyBand.'//trim(ftag)//'.dat')
      ifiei=ifile_handle()
      open(ifiei,file='Band.'//trim(ftag)//'.dat')
      ifidi=ifile_handle()
      open(ifidi,file='BandDiag.'//trim(ftag)//'.dat')
      
!! -----     
      allocate(ovlm(nev,nev),ovlm0(nev,nev),hamm(1:nev,1:nev),fac(napw0,napw0),bqr(napw0,nev),bb(nev),bbx(nev))
      allocate(ovlmi(nev,nev),bqx(napw0),t_zv(1:nev,1:nev),eemat(nev,nev))
      allocate(percent(napw0,ndatx))
      allocate(eee(ndatx,nx),eeex(ndatx))
      print *,'gointo ikp'
      do ikp = 1,ndatx          !-ndatx,ndatx*2
        qp = qpp(:,ikp)
        write(6,"('qp=',3f9.3)") qp
!! bn(r)= \sum_G exp( -i(q+G)r ) exp(- 1/Galpha**2 * (q+G-Gn)**2) =sum_Gi bqr(Gi,iev) *exp(i (q+Gi)r)
        do iev=1,nev
           do i=1,napw0
c             xx(1)=sum(plat(:,1)*(qp+GG(:,i)-ggn(:,iev)))
c             xx(2)=sum(plat(:,2)*(qp+GG(:,i)-ggn(:,iev)))
c             xx(3)=sum(plat(:,3)*(qp+GG(:,i)-ggn(:,iev)))
c             bqr(i,iev) = exp(- sum(xx**2)/Galpha**2 )

              ggabs = sum(GG(:,i)**2)**.5
              alfac= 0.5d0/Grange**2     !good probably best !poor linear depencecy alfac=1 and gbeta=10
              phase=1d0
cccccccccccccccccccccccccccccccc
c              deltar= [0.5,0.5,0.5]
c              pp1= img*2d0*pi* sum((qp+GG(:,i))*plat(:,1))
c              if(abs(pp1)<1d-6) then
c                 phase1=1d0
c              else
c                 phase1= (exp(pp1) - 1d0)/pp1
c              endif
c              pp2= img*2d0*pi* sum((qp+GG(:,i))*plat(:,2))
c              phase2= (exp(pp2) - 1d0)/pp2
c              pp3= img*2d0*pi* sum((qp+GG(:,i))*plat(:,3))
c              phase3= (exp(pp3) - 1d0)/pp3
c              phase = phase1 !*phase2*phase3

c              pp1 = img*2d0*pi* sum((qp+GG(:,i))*plat(:,1)/2d0)
c              phase= 1d0 + exp(pp1)
              
c
c              deltar= [0.5,0.5,0.5]
c              phase= exp( img*2d0*pi* sum((qp+GG(:,i))*deltar) )
c
c$$$c$$$  !! move atom center              
c$$$              if(iev>nev0) then
c$$$                 deltar= [0.25,0.25,0.25]
c$$$                 phase= exp( img*2d0*pi* sum((qp+GG(:,i))*deltar) )
c$$$              elseif(iev>2*nev0) then
c$$$                 deltar= [0.5,0.5,0.5]
c$$$                 phase= exp( img*2d0*pi* sum((qp+GG(:,i))*deltar) )
c$$$              elseif(iev>3*nev0) then
c$$$                 deltar= [0.75,0.75,0.75]
c$$$                 phase= exp( img*2d0*pi* sum((qp+GG(:,i))*deltar) )
c$$$              endif   
cccccc              !alfac= (Galpha + .03*ggabs**2)/(1+0.1*ggabs**2) !0.3 cause strange at 100 for low energy

              phase=1d0
              bqr(i,iev) = exp(- alfac*sum((qp+GG(:,i)-ggn(:,iev))**2)  )*phase

              if(abs(bqr(i,iev)) >1d-2) then
c                write(6,"(' bqrrr ',2i6, ' bqr=',f10.5 )") i, iev, bqr(i,iev)
             endif   
          enddo
c          bqx = bqr(:,iev)
c          write(6,"(' bqrsss ',i6, ' sum max bqr=',2f10.5 )") iev, sum(bqr(:,iev)), maxval(bqr(:,iev))
       enddo
       
!! overlap matrix between CPWs
        ovlm=0d0
        hamm=0d0
        do iev=1,nev
        do jev=1,nev
           do i=1,napw0
              ovlm(iev,jev)= ovlm(iev,jev)+ dconjg(bqr(i,iev))*bqr(i,jev)
              hamm(iev,jev)= hamm(iev,jev)+ (2*pi/alat)**2* sum((qp+GG(:,i))**2)/ 2d0 *dconjg(bqr(i,iev))*bqr(i,jev)
           enddo
        enddo
        enddo
        ovlm0=ovlm
!! Diagonal elements of ovlm is normalized.
        do iev=1,nev
           ooi= ovlm(iev,iev)**.5
        do jev=1,nev
           ooj= ovlm(jev,jev)**.5
           ovlm(iev,jev)= ovlm(iev,jev)/ooi/ooj
           hamm(iev,jev)= hamm(iev,jev)/ooi/ooj
        enddo
        enddo
!!        
        eemat=0d0
        do i=1,nev
           eemat(i,i)=1d0
        enddo
c        print *,'sum1=',sum(ovlm)

!! Diagonalization of ovlm
        ovlmi = ovlm
        call zhev_tk2( nev, ovlmi, eemat, nev , nevout,
     .           evl, t_zv, .false., .false.,.false., ifig) !stock eigenfunctions z
        do i=1,nev
          write(6,"('EigOV=',3f8.3,' ',i5,d13.5)")qp, i,evl(i)
        enddo
        do i=0,nkpe
          if(ikp==ikpe(i)) write(ifp,"('eigovratio=',i5,d13.5)") ikp, evl(1)/evl(nev)
          eigov(i+1)= evl(1)/evl(nev)
        enddo
       write(6,"(2f9.4,i3,  ' ! Grange Gscale cases')") Grange,Gscale, cases


!! Projection weight of GG vectors in the CPW space.
!!   Can we reproduce exp(i q+GG r)?
        ovlmi=ovlm0
        call matcinv(nev,ovlmi)
c       print *,'sum2=',sum(ovlmi),sum(matmul(ovlmi,ovlm))
        gtestmx = gnmax
        write(6,"(' repro:')")
        do i=1,napw0
         GG(:,i) =  matmul(qlat,igv2x0(:,i))
         if(sum((GG(:,i))**2)< gtestmx) then
            bb(:) = bqr(i,:)
            bbx = matmul(ovlmi,bb)
            wsum = sum(bb*bbx)
            percent(i,ikp)=wsum*100
c            write(6,"(' repro: 'i7,x,3f8.3,x,2f8.3,2x,f8.3,' percent')") i,GG(:,i),sum(GG(:,i)**2)**.5,
c     &       (2*pi/alat)**2*sum((qp+GG(:,i))**2)/2d0,  wsum*100
         endif
        enddo

        
!! Empty sphere dispersion
        ix=0
        do i = 1,napw0
           if(sum(GG(:,i)**2)< gnmax*efac) then
              ix=ix+1
              eee(ikp,ix) =(2*pi/alat)**2* sum((qp+GG(:,i))**2)/ 2d0  
c              write(ifieg,"('QPandEG=',f9.4,3f8.3,' ',f10.4)") xqpp(ikp),qp, eee(ikp,ix)
           endif
        enddo
        
!! Diagonalization
        call zhev_tk2( nev, hamm , ovlm , nev , nevout,
     .           evl, t_zv, .false., .false.,.false., ifig) !stock eigenfunctions z
        do i=1,nev
          write(ifiei,"('QPandEigen=',f9.4,3f8.3,' ',f10.4)")xqpp(ikp),qp, evl(i)
        enddo

!! Diagonal only case
        ovlm=0d0
        hamm=0d0
        do iev=1,nev
           jev=iev
           do i=1,napw0
              ovlm(iev,jev)= ovlm(iev,jev)+ bqr(i,iev)*bqr(i,jev)
              hamm(iev,jev)= hamm(iev,jev)+ (2*pi/alat)**2* sum((qp+GG(:,i))**2)/ 2d0 *bqr(i,iev)*bqr(i,jev)
           enddo
           write(ifidi,"('QPandDiag=',f9.4,3f8.3,' ',f10.4)") xqpp(ikp), qp, dreal(hamm(iev,iev)/ovlm(iev,iev))
        enddo
        do ikx = 0,nkpe
           if(ikp==ikpe(ikx)) then
              eeex=eee(ikp,1:nx)
              call bubblesort(nx,eeex)
              diffe=0d0
              do ix = 1,ikk
                 if(ikpix(1,ix)==ikp) then 
                    i = ikpix(2,ix)
c     write(6,*) 'ikpix ',ikpix(1:2,ix)
                    write(ifp,"('Eig=',2i5,x,f9.4,x,f9.4,x,1f10.6)")
     &                   ikp, i, evl(i)*13.605, eeex(i)*13.605, evl(i)*13.605-eeex(i)*13.605
                    if(diffe<abs(evl(i)*13.605-eeex(i)*13.605)) diffe = evl(i)*13.605-eeex(i)*13.605
                 endif
              enddo
              ediff(ikx+1)=diffe
           endif   
        enddo
      enddo

!! Empty sphere dispersion
      do ix=1,nx
      do ikp = 1,ndatx
         write(ifieg,"('QPandEG=',f9.4,3f8.3,' ',f10.4)") xqpp(ikp),qp, eee(ikp,ix)
      enddo
         write(ifieg,*)
      enddo
      close(ifp)

      ifp2=ifile_handle()
      open(ifp2,file='CPWline'//trim(ftag)//'.dat')
      write(ifp2,"(2f7.3,f9.3,i3, d12.3,x,d12.3' ! Grange Gscale gmx cases ovl diffe ')")
     &  Grange,Gscale,gmx, cases, minval(eigov), maxval(abs(ediff)) 
      close(ifp2)
      
      fff1='EmptyBand.'//trim(ftag)//'.dat'
      fff2='Band.'//trim(ftag)//'.dat'
      fff3='BandDiag.'//trim(ftag)//'.dat'

      ifp3=ifile_handle()
      open(ifp3,file='cpwplot.'//trim(ftag)//'.glt')
      aaa="plot '"//trim(fff1)//
     &    "' u ($2):($6*13.605) lt 1 lw 0.5 w l, '"//trim(fff2)//"' u ($2):($6*13.605) lt 2 lw 1"
      write(ifp3,"(a)")trim(aaa)
      write(ifp3,"(a)")"pause -1"
      close(ifp3)

!! Projection weight of GG vectors in the CPW space.
!!   Can we reproduce exp(i q+GG r)?
      ifp4=ifile_handle()
      open(ifp4,file='CPWrepro'//trim(ftag)//'.dat')
      gtestmx = gnmax
      do i=1,napw0
         GG(:,i) =  matmul(qlat,igv2x0(:,i))
         if(sum((GG(:,i))**2)< gtestmx) then
            write(ifp4,*)
            do ikp=1,ndatx
               qp = qpp(:,ikp)
               write(ifp4,"(i7,x,3f8.3,x,2f9.4,2x,f9.4,' percent')") i,GG(:,i),sum(GG(:,i)**2)**.5,
     &              (2*pi/alat)**2*sum((qp+GG(:,i))**2)/2d0, percent(i,ikp)
            enddo
         endif
      enddo
      close(ifp4)
      
      end program lmfCPW

      
      subroutine bubblesort(N,array)
!sikinote, 2016/08/08
      implicit none
      integer,intent(in)::N
      double precision,intent(inout)::array(1:N)
      integer::i,j
      double precision::t
      
      do i=1,N-1
      do j=i+1,N
         if(array(i) .gt. array(j))then
            t=array(i)
            array(i)=array(j)
            array(j)=t
         end if
      end do
      end do

      return
      end subroutine bubblesort
