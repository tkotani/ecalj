      program zeotest 
      real(8) ::ovvpp,emaxx,delt,dq(3)
      integer(4):: nmax,nmx,nn1,nn2,mx,nband,ibd,ifb
      integer(4),allocatable::ilbl(:),iixy(:)
      complex(8),allocatable:: cphi(:,:),cphip(:,:)
      real(8),allocatable:: evl(:),evll(:),qq1(:,:),qq2(:,:)
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

      ibw1=185
      ibw2=210

      open(ifb,file='CHIWV.dat',form='unformatted')
      read(ifb) iline, qx(1:3), qposxx
      read(ifb) iq,nqtot,ib1,ib2,ldim2,ndimh,isp
      write(6,"(100i6)") iq,nqtot,ib1,ib2, ldim2,ndimh,isp
      if(isp/=1) stop 'isp/=1'
      allocate( evl(ib1:ib2), evll(ib1:ib2), 
     &   cphip(1:ldim2,ib1:ib2), cphi(1:ldim2,ib1:ib2)
     &   ,iixy(ib1:ib2),ibidx(ib1:ib2))
      rewind(ifb)
      do ibx=ib1,ib2
        ibidx(ibx)=ibx
      enddo
c      write(ifx2,"(26x, 100i8)")  (ibx, ibx=ibw1,ibw2)
      do 
        read(ifb) iline, qx(1:3), qposxx
        read(ifb) iq,nqtot,ib1,ib2, ldim2,ndimh,isp
        print *
c        print *,' ============= iq=',iq,qq(:,iq)
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
          call connection( cphi,cphip,ldim2, ib1, ib2, iixy) 
          ibidx(ib1:ib2) = iixy(ibidx(ib1:ib2)) !continuous band index
c          ilbl(iixy(ib1:ib2))= ilbl(ib1:ib2)
          write(ifx,"(' Index ! ',100i8)") (ibidx(ibx),ibx=ibw1,ibw2)
        endif

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
      enddo
      end

      subroutine connection( cphi,cphip,ldim2,ib1, ib2, iixy)
      integer(4):: iixy(ib1:ib2),iii,ldim
      complex(8):: cphi(1:ldim2,ib1:ib2),cphip(1:ldim2,ib1:ib2)
      real(8) :: dorm(ib2-ib1+1)
c      print *,' connection'
      ldim=ldim2
      do ibx = ib1,ib2
      do iby = ib1,ib2
       dorm(iby-ib1+1) = sum(abs( 
     &   abs(cphi(1:ldim,iby))**2 - abs(cphip(1:ldim,ibx))**2  )) 
       write(666,*) ibx,iby,dorm(iby-ib1+1)
       enddo
       iii = minloc( dorm,dim=1 )
       iixy(ibx) =  iii + ib1 - 1
      enddo
      end







