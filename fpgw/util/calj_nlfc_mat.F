C -----
C Utinilty routine to calculate Im(chipm)
C set nmbas!
C fname = ChiPM*.nlfc.mat
C -----------------------------
      implicit none
      integer(4),parameter:: nx=1000
      integer(4):: ifx=1011,ifchipm,iw,i,ifmat,j
      real(8):: q(3),jzero,jzero2,  freq(nx), dw=0.2,rrr,iii,polinta !meV
      real(8):: rchipmm(nx),ichipmm(nx),omg,rydberg,mmom,jval,rrrx
     & ,drchi(nx),freqm(nx),czero,omgw,ddd1,ddd2,jdiff,jpp,jpm,eout
     &  ,eout1,etou2,pi,mmomq,eout2,eouts1
     & ,jzero3,polinta1,absq, mmnorm0,bbb,ccc,meff,domega
      character(30) fname,aaa
      character(8) ax,bx
      integer(4)::nmbas,ix,iy
      complex(8),allocatable:: x0mean(:,:,:),x0mat(:,:),x0matinv(:,:)
     &  ,x0matinv0(:,:),xinv(:,:,:)
      real(8),allocatable::
     & r0mean(:,:,:),i0mean(:,:,:),eval(:),r0mean0(:,:),eee(:),
     & freq2(:),eee2(:),a11(:), momsite(:),usite(:),usite2(:),usite3(:),
     & uu0(:),uu1(:),uuo(:,:)
      complex(8):: odd1,odd2,phase,img=(0d0,1d0)
      character(300) :: fff,ggg,hhh
      integer(4),allocatable:: imbas(:)
      real(8),allocatable:: mmnorm(:)

      real(8),allocatable:: usiten(:,:),drchin(:,:)

      complex(8),allocatable:: meffi(:,:), meffix(:,:),dcc(:)

      complex(8),allocatable:: oo(:,:),zz(:,:),jjmat(:,:),sqm(:,:)
      integer(4):: nmx,nev

      real(8),allocatable:: e1(:),e2(:)
      complex(8),allocatable:: eiqrm(:)

      integer(4):: nnmx=60
! This is for ChiPM*.nolfc.mat ! oct2006
      pi=4*atan(1d0)
c      write(6,*)"pi=",pi

c--- Get mmom
c      ifchipm=1015
c      open(ifchipm, file='ChiPM.log')
c      read(ifchipm,*) ax, mmom,bx, mmnorm0
c      write(6,*) ax,mmom
c      if(bx/='svec**2=') stop 'No svec**2 in ChiPM.log'
c      write(6,*) bx,mmnorm0
c      close(ifchipm)

      read(5,'(a)') fff
      ggg= adjustl( fff( 1:index(fff,' '))  )
      hhh= adjustl( fff( index(fff,' '):  )  )
      print *,' ggg=',ggg
      print *,' hhh=',hhh
      read(ggg,*) fname
      ifchipm=1015
      open(ifchipm, file=fname)
      read(ifchipm,*) nmbas
      allocate(imbas(nmbas),momsite(nmbas),mmnorm(nmbas),eiqrm(nmbas),
     &  e1(nmbas),e2(nmbas))
      read(ifchipm,*) imbas(1:nmbas)
      read(ifchipm,*) momsite(1:nmbas)
      read(ifchipm,*) mmnorm(1:nmbas)

      read(ifchipm,*) (e1(ix),e2(ix),ix=1,nmbas)
      eiqrm(1:nmbas) = dcmplx(e1(1:nmbas),e2(1:nmbas))


      if(len(trim(hhh))/=0) then
        print *,' header file=',hhh
        read(hhh,*) fname
        ifchipm=1016
        open(ifchipm, file=fname)
      endif

      write(6,*)'mmnorm=', mmnorm(1:nmbas)

c----------------------------

      mmnorm0= mmnorm(1)**2
      mmom   = momsite(1)

      allocate(x0mean(nmbas,nmbas,0:nx),eee(nx),a11(0:nx),
     &    r0mean(0:nx,nmbas,nmbas) ,r0mean0(nmbas,nmbas),
     & i0mean(0:nx,nmbas,nmbas),x0mat(nmbas,nmbas)
     & ,x0matinv0(nmbas,nmbas),x0matinv(nmbas,nmbas),
     &  usite(nmbas),usite2(nmbas),usite3(nmbas)
     & ,uu0(nmbas),uu1(nmbas),xinv(nmbas,nmbas,nx))
       allocate(usiten(nx,nmbas),drchin(nx,nmbas)
     & ,uuo(nmbas,nx)  )
       uu1=0d0; uuo=0d0
      allocate(meffi(nmbas,nmbas),dcc(nx))
      allocate(meffix(nmbas,nmbas))

cccccccccccccccccccccccccccccc
C this is the case when svec is used.
c      momsite=1d0
c      momsite(3)=0d0
c      momsite(4)=0d0
cccccccccccccccccccccccccccccc

c      read(ifchipm,*)q, freq(1)
c      if(sum(abs(q(:)))==0d0.and.freq(1)==0d0) then
c        ifx = 1010 
c      endif
c      rewind ifchipm
c      read(ifchipm,*) (momsite(ix),ix=1,nmbas)

! *.nolfc.mat---> Now not 1st line skip
!      read(ifchipm,*)
      do iw=1,nnmx
        read(ifchipm,*)q, freq(iw)
     &   , ( (r0mean(iw,ix,iy),i0mean(iw,ix,iy),ix=1,nmbas),iy=1,nmbas)
        x0mean(:,:,iw)= dcmplx(r0mean(iw,:,:),i0mean(iw,:,:))

        if(iw<4.or.mod(iw,10)==1)  then
         write(6,'("read ",3f8.4,2x,f13.7,2x,255d12.4)') 
     &   q, freq(iw) !,x0mean(iw,1:nmbas,1:nmbas)
         do iy=1,nmbas
          write(6,'($,"      ",24x,2x,255d12.4)') 
     &    (dreal(x0mean(ix,iy,iw)),ix=1,nmbas)
          write(6,'(" | ",255d12.4)') 
     &    (dimag(x0mean(ix,iy,iw)) ,ix=1,nmbas)
         enddo
        endif
C--- x0mean is now normaized as
C    x0mean = <m|x0|m'>/sqrt(<m|m>) /sqrt(<m'|m'>)
        print *,' iw x0mean=',iw,x0mean(:,:,iw)
        call matcinv(nmbas, x0mean(:,:,iw))

C ... determine usite and uu0 matrix.
        if(sum(abs(q(:)))==0d0) then
         if(iw==1) then
          ifx=1010 
          usite(1:nmbas)=
     &    matmul(dreal(x0mean(1:nmbas,1:nmbas,iw)),mmnorm(1:nmbas))
          write(6,'(" usite = ",255d12.4)')  usite(1:nmbas)

          uu0=0d0
          do i=1,nmbas
            if(abs(momsite(i))>1d-3 ) uu0(i)=usite(i)/mmnorm(i)
          enddo
          write(6,*)
c          do iy=1,nmbas
c           write(6,'("      ",24x,2x,255d12.4)') 
c     &     (dreal(x0mean(ix,iy,iw)) - uu0(ix,iy),ix=1,nmbas)
c          enddo
c         write(6,'(" usite= ",255d12.4)')  (usite(ix) ,ix=1,nmbas)
          write(6,'(" uu0 = ",255d12.4)') 
     &     (uu0(ix) ,ix=1,nmbas)
         endif 

C ... determine uu1 matrix.
ccc         if(iw<=10) then
          usiten(iw,1:nmbas) = 
     &    matmul(dreal(x0mean(1:nmbas,1:nmbas,iw)),mmnorm(1:nmbas))
c          write(6,"(' usiten =',12f23.16)")  usiten(iw,:)
c          do i=1,nmbas
c          usiten(iw,i) = dreal(x0mean(i,i,iw))*mmnorm(i)
c          enddo
c         do iy=1,nmbas
c          write(6,'("vvv      ",24x,2x,255d12.4)') 
c     &    (dreal(x0mean(ix,iy,iw)),ix=1,nmbas)
ccc         endif

         if(iw==10) then
          do ix=1,9
            drchin(ix,:) 
     &       = (usiten(ix+1,:) - usiten(ix,:))/(freq(ix+1)-freq(ix))  !Feb24-2007
            freqm(ix) = 0.5d0*(freq(ix+1)+freq(ix))
            write(6,"(' :f drchi=',12f23.16)")  freqm(ix),drchin(ix,:)
          enddo
          do i=1,nmbas
            omg  = 0d0 !freqm(8) !0d0
            rrr  = polinta(omg, freqm(2:8),drchin(2:8,i),7)
ccccccccccccccccccccc
            print *,'kkk : rrr=',rrr
ccccccccccccccccccccc
            if(abs(momsite(i))>1d-3 ) then
              uu1(i)= ( - momsite(i)/mmnorm(i) + rrr) /mmnorm(i)
           endif
          enddo
          omg =   freq(iw)     - freq(1)
          write(6,'(" kkk: uu1= ",i3,255d12.4)')iw,omg,(omg*uu1(ix), ix=1,nmbas)
          write(6,'(" mom/mmnorm**2 = ",255d12.4)') 
     &     (momsite(ix)/mmnorm(ix)**2,ix=1,nmbas)
         endif 

c         if(iw==10) then
c          ix= iw
c          do ix=9,9
c            drchin(ix,:) 
c     &       = (usiten(9,:) - usiten(1,:)) !/(freq(ix+1)-freq(ix))  !Feb24-2007
c            freqm(ix) = 0.5d0*(freq(ix+1)+freq(ix))
c            write(6,"(' :f drchi=',12f23.16)")  freqm(ix),drchin(ix,:)
c          enddo
          omg =   freq(iw)     - freq(1)
          do i=1,nmbas
            rrr =  (usiten(iw,i) - usiten(1,i)) 
            if(abs(momsite(i))>1d-3 ) then
              uuo(i,iw)= ( - omg*momsite(i)/mmnorm(i) + rrr) /mmnorm(i)
           endif
          enddo
          write(6,'(" kkk: uuo= ",i3,255d12.4)') iw,omg,(uuo(ix,iw), ix=1,nmbas)
          write(6,'(" kk : mom/mmnorm**2 = ",255d12.4)') 
     &     (momsite(ix)/mmnorm(ix)**2,ix=1,nmbas)
c         endif 

c         if(iw==2) then
c          usite2 = 
c     &    matmul(dreal(x0mean(1:nmbas,1:nmbas,iw)),mmnorm(1:nmbas))
c         endif
c         if(iw==3) then
c          usite3 = 
c     &    matmul(dreal(x0mean(1:nmbas,1:nmbas,iw)),mmnorm(1:nmbas))
c          uu1=0d0
c          domega = (freq(3)-freq(2))
c          do i=1,nmbas
c            if(abs(momsite(i))>1d-3 ) then
c              uu1(i,i)= ( - momsite(i)/mmnorm(i) 
c     &                    + (usite3(i)-usite2(i)) /domega ) /mmnorm(i)
c           endif
c          enddo
c          write(6,'(" uu1 = ",255d12.4)') (uu1(ix,ix), ix=1,nmbas)
c          write(6,'(" mom/mmnorm**2 = ",255d12.4)') 
c     &     (momsite(ix)/mmnorm(ix)**2,ix=1,nmbas)
c         endif 
        endif

      enddo
C .. write or read uu0uu1
      open(ifx,file='uu0uu1')
      if(ifx==1010) then
        write(ifx,'(255d23.15)') (uu0(ix),ix=1,nmbas)
        write(ifx,'(255d23.15)') (uu1(ix),ix=1,nmbas)
        do iw=1,nnmx
        write(ifx,'(255d23.15)') (uuo(ix,iw),ix=1,nmbas)
        enddo
      elseif(ifx==1011) then
         print *,' read uu0 uu1'
        uu0=0d0; uu1=0d0
        read(ifx,*) (uu0(ix),ix=1,nmbas)
        read(ifx,*) (uu1(ix),ix=1,nmbas)
        do iw=1,nnmx
          read(ifx,'(255d23.15)') (uuo(ix,iw),ix=1,nmbas)
        enddo
      endif
      close(ifx)
      print *,' use uu0 =',(uu0(ix),ix=1,nmbas)
      print *,' use uu1 =',(uu1(ix),ix=1,nmbas)

C --- get xinv = x
      allocate(eval(nmbas))
      do iw=1,nnmx
cccccccccccccccccccccccccc
c        xinv(:,:,iw) = x0mean(:,:,iw) - uu0(1:nmbas,1:nmbas) 
c     &   - freq(iw) * uu1(1:nmbas,1:nmbas) 
        xinv(:,:,iw) = x0mean(:,:,iw) 
        do ix=1,nmbas
          xinv(ix,ix,iw) = xinv(ix,ix,iw) -  uu0(ix)-uuo(ix,iw) 
        enddo
cccccccccccccccccccccc

        xinv(:,:,iw) = xinv(:,:,iw) + transpose(dconjg(xinv(:,:,iw)))
        xinv(:,:,iw) = .5d0*xinv(:,:,iw)
cccccccccccccccccccccccccccccc
         write(6,'($,"xinv ",3f8.4,2x,f13.7,2x,255d12.4)') 
     &   q, freq(iw) !,x0mean(iw,1:nmbas,1:nmbas)
         do iy=1,nmbas
          write(6,'($,255d12.4)') 
     &    (dreal(xinv(ix,iy,iw)),ix=1,nmbas)
c          write(6,'($," | ",255d12.4)') 
          write(6,'($,255d12.4)') 
     &    (dimag(xinv(ix,iy,iw)),ix=1,nmbas)
         enddo
ccccccccccccccccccccccccccccccccccccc
        call diagno00(nmbas,xinv(:,:,iw),eval)
         write(6,*)
c        write(6,"(' eval  ', 3f8.4,2x,i4,2x,f13.7,2x,255f14.8)") 
c     &  q, iw, rydberg()*1d3*freq(iw), (eval(ix),ix=1,nmbas)
        eee(iw) = maxval(eval) ! &   , ( (x0mat(ix,iy),ix=1,nmbas),iy=1,nmbas)
        if(iw==1) then
          jpp=  x0mean(1,1,iw)/4d0
          jpm= abs(-x0mean(1,2,iw)/4d0)
          write(6,*) 'jpp=',jpp
          write(6,*) 'jpm=',jpm
          x0matinv0= xinv(:,:,iw)
          jzero=1d99
          if(sum(abs(q(:)))==0d0) jzero = jpp + abs(jpm)
        endif
        a11(iw) = xinv(1,1,iw)
      enddo

cccccccccccccccccccccccccccccccccccccccccccccc
c      do i=1,nmbas
c      do j=1,nmbas
c      if(i/=j) meffi(i,j)=0d0
c      enddo
c      enddo
cccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C 2x2 matrix case
c      close(ifmat)
! effective  1/m = mmnorm0/mmom**2 * rrr
! mmom**2/(rrr*mmnorm0)
! Ferro case
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 1092 format(f19.15,2x,3f19.15, f10.2,a)
 1093 format(f19.15,2x,3f19.15, f10.2,d13.5,a)
c-----------
      absq=sqrt(sum(q**2))

C ... for pole F case
c      if(eee(1) <1d-10) then
        eout=0d0
        rrrx=-1d0
        do 
        eout = eout + 1d-6/rydberg()
        rrrx  = polinta(eout, freq, eee, nnmx)
        if(rrrx>0) exit
        enddo
        eouts1=eout
c      else
c       eout= 4d0/rrr * (-jpp+jzero-abs(jpm)) 
c       write(6,1092) absq, q,  mmom,
c     &  rydberg()*1d3*eout 
c      endif


C... determine meffi
       omg= eouts1
c      omg=40d0/ rydberg()/1d3
c      omg   = 0d0
      do i=1,nmbas
      do j=1,nmbas
        do iw=1,9
        dcc(iw)  = (xinv(i,j,iw+1) - xinv(i,j,1))/(freq(iw+1)-freq(1))
        freqm(iw) = 0.5d0*(freq(iw+1)+freq(iw))
        enddo
        meffi(i,j)  = polinta(omg, freqm(1:8),dreal(dcc(1:8)),8)
     &          +img* polinta(omg, freqm(1:8),dimag(dcc(1:8)),8)
      enddo
      enddo
      do iw=1,9
        drchi(iw) = (a11(iw+1) - a11(iw))/(freq(iw+1)-freq(iw))  !Feb24-2007
        freqm(iw) = 0.5d0*(freq(iw+1)+freq(iw))
        write(6,"(' f drchi=',2f23.16)")  freqm(iw),drchi(iw)
      enddo
      rrr  = polinta(omg, freqm(2:8),drchi(2:8),7)
      write(6,"(2f23.16,' ! drchi=dchipm^-1/domeg \approx 1/m')") 
     & rrr, 1d0/rrr
!
      do i=1,nmbas
      do j=1,nmbas
        x0matinv0(i,j) =  x0matinv0(i,j)
     &       * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
        meffi(i,j)     =  meffi(i,j) 
     &       * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
      enddo
      enddo
      meff = mmom**2/(rrr*mmnorm0)
      print *,' meff=',meff




C ... for pole search AF case
      allocate(freq2(-100:100),eee2(-100:100))
      eee2(1:100)=eee(1:100)
      freq2(1:100)=freq(1:100)
      do i=2,nnmx
        eee2(2-i) =  eee(i)
        freq2(2-i)= -freq(i)
      enddo
      do iw=-10,10
        write(6,"('eval2  ',i3,2f12.6)")
     &   iw,rydberg()*1d3*freq2(iw), eee2(iw)
      enddo
      eout=0d0
      rrrx=-1d0
      do 
        eout = eout + 1d-7/rydberg()
        rrrx  = polinta(eout, freq2(-28:), eee2(-28:), 59)
        if(rrrx>0) exit
      enddo
      write(6,*)
      write(6,*) '========= RESULTS SW (mRy) ============'
      write(6,*)" --- AF case ---"
      write(6,1092) 
     &  absq, q, rydberg()*1d3*eout , 
     &  ' !   pole search2 quadratic AF case'
      deallocate(freq2,eee2)
 1010 continue
c      write(6,1093) absq,q, mmom,
c     &  rydberg()*1d3* 4d0/rrr * sqrt( (jpp-jzero)**2-jpm**2)
c     &  ,(jpp-jzero)**2-jpm**2
c     & ,' !   rrr slope AF'
!AF case
      eout2 = meff*sqrt(x0matinv0(1,1)*x0matinv0(2,2)- x0matinv0(1,2)*x0matinv0(2,1))
      write(6,1092)
     &  absq, q,
     &  rydberg()*1d3* eout2,
     &  ' ! diagonalization Formula AF 2atom'
c      jdiff= ( (jpp - jzero)**2 - jpm**2 )* mmnorm0**2
c      omgw = rydberg()*1d3* 4d0/mmom* sqrt(jdiff)
c      write(6,1092)
c     &  sqrt(sum(q**2)),q,  omgw,
c     &  ' ! jpp-jpm formula AF 2atom'
c      write(6,*)

C--------------------------------------------
      write(6,*)" --- Ferro case ---"
        write(6,1092)
     &  absq, q, rydberg()*1d3*eouts1,
     &  ' !q  m  omega(meV). pole search1'
c      write(6,1093) absq,q, mmom,
c     &  rydberg()*1d3* 4d0/rrr * (-jpp+jzero-abs(jpm))
c     &  ,(jpp-jzero) !**2-jpm**2
c     & ,' !   rrr slope Ferro'
      bbb = x0matinv0(1,1) + x0matinv0(2,2)
      ccc = x0matinv0(1,1)*x0matinv0(2,2)- x0matinv0(1,2)*x0matinv0(2,1) 
      eout1 = meff*(-bbb - sqrt(bbb**2-4*ccc))/2d0
      write(6,1092)
     &  absq, q, 
     &  rydberg()*1d3* eout1, ' ! diagonalization Formula F 2atom'

c      jdiff= (-jpp+jzero-abs(jpm))
c      omgw = rydberg()*1d3* 4d0/mmom* jdiff
c      write(6,1092) absq,q, omgw,
c     &  ' ! jpp-jpm formula F 2atom'
      mmomq=mmom**2/(rrr*mmnorm0)
      write(6,*)
      write(6,"(f14.7,' !!static  moment M(0) ')") mmom
      write(6,"(f14.7,' !!               M(q) ')") mmomq
      write(6,"(f14.7,' !! No uu1 case        ')") 
     &   mmom**2/((rrr+uu1(1))*mmnorm0)
c1d0/rrr, rrr*mmnorm0
      write(6,"(f14.7,' !!Stoner I(eV)     ',)") 
     &     rydberg()*4*jzero *mmnorm0 /mmom**2
      write(6,*)
      ifmat= 6!1017
c      open(ifmat, file="Jmat")!, access="append")
      write(ifmat,"('JMAT: ',3d18.10, 3x, 255d18.10)") 
     &  q, ((x0matinv0(ix,iy),ix=1,nmbas),iy=1,nmbas) 
      write(ifmat,"('MMAT: ',3d18.10, 3x, 255d18.10)") 
     &  q, ((meffi(ix,iy),ix=1,nmbas),iy=1,nmbas) 

C---------------
      write(ifmat,"('Mmm: ',3d18.10, 3x, 255d18.10)") 
     &  q, ((1d0/meffi(ix,iy),ix=1,nmbas),iy=1,nmbas) 

c        phase = exp( 2*pi*img*sum(q*(/0d0,-0.5d0,-0.5d0/)) )
c        write(6,'("JMM2 ",3f8.4,x,f13.7,x, 4(" (",2d12.4,") ") ') 
c     &   q, freq(iw) 
c     &   ,x0matinv0(1,2) !*dconjg(phase)
c     &   ,x0matinv0(1,2)*phase
c     &   ,x0matinv0(1,2)*dconjg(phase)
c      write(6,'(" phase",2f18.4,"  ",2d13.5,3f8.4)') 
c     &  imag( log(x0matinv0(1,2)*phase) )/(3.1415926),
c     &  2*sum(q*(/0d0,-0.5d0,-0.5d0/)),
c     &  (x0matinv0(1,2)*phase), q

c     &  sum(q*(/0d0,-0.5d0,-0.5d0/)), 
c     &  (x0matinv0(1,2)),q
c     &  dimag(exp(2*pi*img*sum(q*(/1d0,1d0,1d0/))))/
c     &  dreal(exp(2*pi*img*sum(q*(/1d0,1d0,1d0/)))),
c     &  dimag(x0mean(1,2,iw))/dreal(x0mean(1,2,iw))

      allocate(   oo(nmbas,nmbas),zz(nmbas,nmbas),
     &         jjmat(nmbas,nmbas),sqm(nmbas,nmbas))
       oo=0d0
       do i=1,nmbas
        oo(i,i)=1d0
       enddo
       nmx=nmbas
       nev=nmbas
ccccccccccccccccccccc
c       do ix=1,nmbas
c       do iy=1,nmbas
c         if(ix/=iy) meffi(ix,iy)=0d0
c       enddo
c       enddo
ccccccccccccccccccc
       meffix=meffi
       call diagcv(oo,meffix,zz,nmbas, eval,nmx,1d99,nev)
       write(6,*)' oo eig=',eval
       oo=0d0
       do i=1,nmbas
        oo(i,i)= 1d0/sqrt(eval(i))
       enddo
       sqm   = matmul(zz, matmul(oo, transpose(dconjg(zz))) )
       jjmat = matmul(sqm, matmul(x0matinv0, sqm))
       do ix=1,nmbas
       do iy=1,nmbas
       jjmat(ix,iy)= jjmat(ix,iy)/sqrt(momsite(ix)*momsite(iy))
       enddo
       enddo
      write(ifmat,"('JJMAT: ',3d18.10, 3x, 255d18.10)") 
     &  q, (( jjmat(ix,iy) ,ix=1,nmbas),iy=1,nmbas) 

      call diagno00(nmbas,jjmat,eval)
      print *,' momsite=',momsite
      write(6,"('e sw  ', 3f8.4,2x,255f14.8)") 
     &   q,(-momsite(1)*1d3*rydberg()*eval(ix),ix=1,nmbas)
C----------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccc
       nmx=nmbas
       nev=nmbas
       meffix=meffi
       call diagcv(meffix,x0matinv0,zz,nmbas, eval,nmx,1d99,nev)
      write(6,"('t sw  ', 3f8.4,2x,255f14.8)") 
     &  q,(-1d3*rydberg()*eval(ix),ix=1,nmbas)
      write(6,"('e mom ', 3f8.4,2x,255f14.8)") 
     &  q,(sum( dconjg(zz(:,ix))*zz(:,ix)),ix=1,nmbas)
ccccccccccccccccccccccccccccccccccccccccc


      end








      real(8) function rydberg()
      rydberg=13.6058d0
      end

      double precision function polinta1(x,xa,ya,n)
c----------------------------------------------------------------------
c     Given arrays xa and ya, each of length n and given value x,
c     this function returns a value polint. If p(x) is the polynominal
c     of degree ndg such that p(xa(i))=ya(i), i=ns,..,ns+ndg then
c     the returned value polint=p(x). ns is obtained by hunting.
c     See Numerical Recipes
c     coded by H.Akai
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (ndgmx=1, nmx=ndgmx+1)
      dimension xa(n),ya(n),c(nmx),d(nmx)
      logical ascnd
      save jlo
      data jlo/0/ , small/1d-30/
      ndg=min(ndgmx,n-1)
      ndt=ndg+1
      ascnd=xa(n) .gt. xa(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
      jlo=0
      jhi=n+1
      go to 30
      endif
      inc=1
      if(x .gt. xa(jlo) .eqv. ascnd) then
   10 jhi=jlo+inc
      if(jhi .gt. n) then
      jhi=n+1
      else if(x. gt. xa(jhi) .eqv. ascnd) then
      jlo=jhi
      inc=inc+inc
      go to 10
      endif
      else
      jhi=jlo
   20 jlo=jhi-inc
      if(jlo .lt. 1) then
      jlo=0
      else if(x .lt. xa(jlo) .eqv. ascnd) then
      jhi=jlo
      inc=inc+inc
      go to 20
      endif
      endif
   30 if(jhi-jlo .ne. 1) then
      jm=(jhi+jlo)/2
      if(x .gt. xa(jm) .eqv. ascnd) then
      jlo=jm
      else
      jhi=jm
      endif
      go to 30
      endif
      nlo=max(1,jlo-ndg/2)
      nhi=min(n,nlo+ndg)
      nlo=nhi-ndg
      if(jlo .eq. 0) then
      ns=1
      else if(jlo .eq. n) then
      ns=ndt
      else if(abs(x-xa(jlo)) .lt. abs(x-xa(jhi))) then
      ns=jlo-nlo+1
      else
      ns=jhi-nlo+1
      endif
      do 40 i=1,ndt
      ii=nlo+i-1
      c(i)=ya(ii)
   40 d(i)=ya(ii)
      polint=ya(nlo+ns-1)
      ns=ns-1
      do 60 m=1,ndg
      do 50 i=1,ndt-m
      ii=nlo+i-1
      ho=xa(ii)-x
      hp=xa(ii+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
c
c     an error can occur if two xa's are identical
      if(abs(den) .lt. small) then
      write(6,1000)
 1000 format('   ***wrn in polint...data error')
      stop
      endif
c
      den=w/den
      d(i)=hp*den
   50 c(i)=ho*den
      if(2*ns .lt. ndt-m) then
      dy=c(ns+1)
      else
      dy=d(ns)
      ns=ns-1
      endif
   60 polint=polint+dy
c takao
      polinta1=polint
      return
      end



      double precision function polinta(x,xa,ya,n)
c----------------------------------------------------------------------
c     Given arrays xa and ya, each of length n and given value x,
c     this function returns a value polint. If p(x) is the polynominal
c     of degree ndg such that p(xa(i))=ya(i), i=ns,..,ns+ndg then
c     the returned value polint=p(x). ns is obtained by hunting.
c     See Numerical Recipes
c     coded by H.Akai
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (ndgmx=4, nmx=ndgmx+1)
      dimension xa(n),ya(n),c(nmx),d(nmx)
      logical ascnd
      save jlo
      data jlo/0/ , small/1d-30/
      ndg=min(ndgmx,n-1)
      ndt=ndg+1
      ascnd=xa(n) .gt. xa(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
      jlo=0
      jhi=n+1
      go to 30
      endif
      inc=1
      if(x .gt. xa(jlo) .eqv. ascnd) then
   10 jhi=jlo+inc
      if(jhi .gt. n) then
      jhi=n+1
      else if(x. gt. xa(jhi) .eqv. ascnd) then
      jlo=jhi
      inc=inc+inc
      go to 10
      endif
      else
      jhi=jlo
   20 jlo=jhi-inc
      if(jlo .lt. 1) then
      jlo=0
      else if(x .lt. xa(jlo) .eqv. ascnd) then
      jhi=jlo
      inc=inc+inc
      go to 20
      endif
      endif
   30 if(jhi-jlo .ne. 1) then
      jm=(jhi+jlo)/2
      if(x .gt. xa(jm) .eqv. ascnd) then
      jlo=jm
      else
      jhi=jm
      endif
      go to 30
      endif
      nlo=max(1,jlo-ndg/2)
      nhi=min(n,nlo+ndg)
      nlo=nhi-ndg
      if(jlo .eq. 0) then
      ns=1
      else if(jlo .eq. n) then
      ns=ndt
      else if(abs(x-xa(jlo)) .lt. abs(x-xa(jhi))) then
      ns=jlo-nlo+1
      else
      ns=jhi-nlo+1
      endif
      do 40 i=1,ndt
      ii=nlo+i-1
      c(i)=ya(ii)
   40 d(i)=ya(ii)
      polint=ya(nlo+ns-1)
      ns=ns-1
      do 60 m=1,ndg
      do 50 i=1,ndt-m
      ii=nlo+i-1
      ho=xa(ii)-x
      hp=xa(ii+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
c
c     an error can occur if two xa's are identical
      if(abs(den) .lt. small) then
      write(6,1000)
 1000 format('   ***wrn in polint...data error')
      stop
      endif
c
      den=w/den
      d(i)=hp*den
   50 c(i)=ho*den
      if(2*ns .lt. ndt-m) then
      dy=c(ns+1)
      else
      dy=d(ns)
      ns=ns-1
      endif
   60 polint=polint+dy
c takao
      polinta=polint
      return
      end


c---------------------------------------------------------------------
      subroutine diagno00(nbloch,wpvcin,eval)
      implicit none
      integer(4):: nbloch,nmx,nev,i
      complex(8),allocatable:: ovlpc(:,:),evecc(:,:)
      real(8)::emx
      complex(8)::wpvc(nbloch,nbloch),wpvcin(nbloch,nbloch)
      real(8):: eval(nbloch)
      allocate( ovlpc(nbloch,nbloch),
     &          evecc(nbloch,nbloch))
      ovlpc = 0d0
      wpvc=wpvcin
      do i = 1,nbloch
        ovlpc(i,i)= 1d0
      enddo
      eval  = 0d0
      nev  = nbloch
      nmx  = nbloch
c      print *,' goto diagcv-----------'
      call diagcv(ovlpc,wpvc, evecc, nbloch, eval, nmx, 1d99, nev)
      deallocate(ovlpc,evecc)
      end
