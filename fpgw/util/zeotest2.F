      program zeotest2 
      real(8) ::ovvpp,emaxx,delt,dq(3)
      integer(4):: nmax,nmx,nn1,nn2,mx,nband,ibd,ifb
      integer(4),allocatable::ilbl(:),iixy(:)
      complex(8),allocatable:: cphi(:,:),cphip(:,:)
      real(8),allocatable:: evl(:),evll(:),qq1(:,:),qq2(:,:)
     &,eig(:),eigp(:)
      integer(4)::ib1,ib2
      integer(4),allocatable:: nqq(:),ibidx(:)
      real(8)::qq(3,1000),q_p(3),qpos(0:1000) ,qx(3)
c
C   ... Read SYML file
        nlinemax = 50
        allocate(nqq(nlinemax),qq1(1:3,nlinemax),qq2(1:3,nlinemax))
        ifsyml=101
        open(ifsyml,file='SYML')
        nline = 0
        do
          nline = nline + 1
          read(ifsyml,*,err=601,end=601)
     &      nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
        enddo
  601   continue
        nline = nline - 1
        iline = 1
        ix=0
        qpos(0)=0d0
        do iline=1,nline
          dq = (qq2(1:3,iline)-qq1(1:3,iline) )/(nqq(iline)-1)
          do ii= 0,nqq(iline)-1
            ix = ix+1 
            qq(1:3,ix) = qq1(1:3,iline) + dq*ii 
            if (ix == 1) q_p=qq(:,ix)
            qpos(ix) = qpos(ix-1) + sqrt( sum((qq(:,ix)-q_p)**2) )
            q_p = qq(:,ix)
          enddo
        enddo
c----------------------------------------------
      ib1 = 185
      ib2 = 210
      ifb = 5000
      ifx =6
      ifx2=0
      ifx3=3000

      ibw1=193
      ibw2=210

      open(ifb,file='CHIWV.dat',form='unformatted')
c      read(ifb) iline, qx(1:3), qposxx
      read(ifb) iq,nqtot,ib1,ib2,ldim2,ndimh,isp
      write(6,"(100i6)") iq,nqtot,ib1,ib2, ldim2,ndimh,isp
      if(isp/=1) stop 'isp/=1'
      allocate( evl(ib1:ib2), evll(ib1:ib2), 
     &   cphip(1:ldim2,ib1:ib2), cphi(1:ldim2,ib1:ib2)
     &   ,eig(ib1:ib2), eigp(ib1:ib2)
     &   ,iixy(ib1:ib2),ibidx(ib1:ib2))
      rewind(ifb)
      do ibx=ib1,ib2
        ibidx(ibx)=ibx
      enddo
c      write(ifx2,"(26x, 100i8)")  (ibx, ibx=ibw1,ibw2)
      do 
c        read(ifb) iline, qx(1:3), qposxx
        read(ifb) iq,nqtot,ib1,ib2, ldim2,ndimh,isp
        print *
c        print *,' ============= iq=',iq,qq(:,iq)
        qx=qq(:,iq)
        print *,' ============= iq=',iq,qx(:)
        read(ifb) evl(ib1:ib2)
        read(ifb) cphi(1:ldim2,ib1:ib2)
c
        write(ifx,"(' Index ! ',100i8)")   (ibx,     ibx=ibw1,ibw2)
        write(ifx,"(' Energy! ',100f8.4)") (evl(ibx),ibx=ibw1,ibw2)
        write(ifx,"(' Wgtsm ! ',100f8.4)") 
     &    (sum(abs(cphi(1:ldim2,ibx))**2),ibx=ibw1,ibw2)

c        do ldim=1,ldim2
c        write(ifx,"(100f8.4)") (abs(cphi(ldim,ibx)),ibx=ibw1,ibw2)
c        enddo

        if(iq>=2) then
          call connection( iq,cphi,cphip,ldim2, ib1, ib2, evl,eigp,
     & iixy) 
          ibidx(ib1:ib2) = iixy(ibidx(ib1:ib2)) !continuous band index
c          ilbl(iixy(ib1:ib2))= ilbl(ib1:ib2)
          write(ifx,"(' Indexm! ',100i8)") (iixy(ibx),ibx=ibw1,ibw2)
          write(ifx,"(' Index ! ',100i8)") (ibidx(ibx),ibx=ibw1,ibw2)
        endif
        write(ifx,"(' Indexsum ! ',100i8)") sum(ibidx(ib1:ib2))
     &  ,sum(iixy(ib1:ib2))

c        ibxx(ib1:ib2)=ilbl(ib1:ib2)
c        evll(ilbl(ib1:ib2)) = evl(ib1:ib2) ! ilbl(ibx)

c-------------------------------
        do ibr=ibw1,ibw2
          write(ifx2+ibr,"(2x,f9.5,3f9.5,2x,f10.5,i4)")
     &    qpos(iq),qx(1:3),evl(ibidx(ibr)),ibidx(ibr)
          write(ifx3+ibr,"(2x,f9.5,3f9.5,2x,f10.5,i4)")
     &    qpos(iq),qx(1:3),evl(ibr)
        enddo
c
        if(iq==nqtot) exit
        cphip = cphi
        eigp = evl
      enddo
      end

      subroutine connection( iq,cphi,cphip,ldim2,ib1, ib2,eig,eigp, 
     &  iixy)
      integer(4):: iixy(ib1:ib2),iii,ldim,iq
      complex(8):: cphi(1:ldim2,ib1:ib2),cphip(1:ldim2,ib1:ib2),
     &  phase
      real(8) :: dorm(ib2-ib1+1),eig(ib1:ib2),eigp(ib1:ib2)
      integer(4)::ifil(ib2-ib1+1)

c      print *,' connection'
      ldim=ldim2

ccccccccccccccccccccccc
c      ibx=185
c      iby=188
c      if(iq>=34.and.iq<=35) then
c        do idim=1,ldim
c         if( abs(cphi(idim,ibx))>0.1.or. abs(cphi(idim,iby))>0.1) then
c         write(8000+iq,"(i4,2f10.3,3x,f10.3)") 
c     &    idim,cphi(idim,ibx),abs(cphi(idim,ibx))
c         write(8000+iq,"(i4,2f10.3,3x,f10.3)") 
c     &    idim,cphi(idim,iby),abs(cphi(idim,iby))
c         write(8000+iq,*)
c         endif 
c        enddo
c      endif
cccccccccccccccccccccc

      ifil=0
      do ibx = ib1,ib2
         idimx = maxloc(abs(cphip(:,ibx)),dim=1)
      do iby = ib1,ib2
c       dorm(iby-ib1+1) = 0d0
cc       do idim=1,ldim
c         dorm(iby-ib1+1) = dorm(iby-ib1+1) + 
c     &     (abs(cphi(1:ldim,iby)) - abs(cphip(1:ldim,ibx)))**2  
c       enddo
c
c         phase=0d0
c         do idim=1,ldim
c          if(abs(cphip(idim,ibx))>0.01d0) then
c           phase = phase + dconjg(cphip(idim,ibx))*cphi(idim,iby)
c          endif 
c         enddo

         phase = sum(dconjg(cphip(1:ldim,ibx))*cphi(1:ldim,iby))
         phase = phase/abs(phase)

         dorm(iby-ib1+1) = 0d0
cccccccccccccccccccccccc
c          if((iq==7.or.iq==8).and.ibx==210.and.iby==210) then
c           write(667,"(' iq=',i4,3d11.3)") iq,phase,abs(phase)
c          endif
cccccccccccccccccccccccc
         do idim=1,ldim
          if(abs(cphip(idim,ibx))>0.01d0) then
            dorm(iby-ib1+1) = dorm(iby-ib1+1) + 
     &      abs(cphi(idim,iby) - cphip(idim,ibx)*phase)**2
          endif
          if(abs(eig(iby)-eigp(ibx))>0.01d0) 
     &    dorm(iby-ib1+1) = 1d20
cccccccccccccccccccccccc
c          if((iq==7.or.iq==8).and.ibx==210.and.iby==210) then
c           write(667,"(' idim=',i4,2x,3d11.3)") 
c     &     idim, cphi(idim,iby),  abs(cphi(idim,iby))
c           write(667,"(' idim=',i4,2x,3d11.3)") 
c     &     idim, cphip(idim,ibx), abs(cphip(idim,ibx))
c           write(667,"(' idim=',i4,2x,3d11.3)") 
c     &     idim, cphip(idim,ibx)*phase,
c     &     abs(cphip(idim,ibx)*phase)
c           write(667,*)
c          endif
ccccccccccccccccccccccccc
         enddo

c         dorm(iby-ib1+1) = sum(abs(
c    &        abs(cphi (1:ldim2,iby))
c     &       -abs(cphip(1:ldim2,ibx))  )**2)
       write(666,"('iq=',i3,2x,2i4,2x,d11.3)") 
     & iq, ibx,iby,dorm(iby-ib1+1)
       enddo

 200   continue
       iii = minloc( dorm,dim=1 )
       if(ifil(iii)==0) then
         iixy(ibx) =  iii + ib1 - 1
         ifil(iii) = 1
       else
         dorm(iii)=1d10
         goto 200
       endif

      enddo

      end











