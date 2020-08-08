      implicit none
      complex(8) :: x0matinv0(2,2)
      real(8):: eout1,bbb,ccc,pi,meff,rydberg

      x0matinv0(1,1)=  -0.1822466259D-02  
      x0matinv0(2,1)=   0.1662788373D-02 
      x0matinv0(1,2)=   0.1662788373D-02  
      x0matinv0(2,2)=  -0.1819846748D-02

      meff= 3.0840716 

      bbb  = x0matinv0(1,1) + x0matinv0(2,2)
      ccc  = x0matinv0(1,1)*x0matinv0(2,2)- x0matinv0(1,2)*x0matinv0(2,1) 
      eout1 = meff*(-bbb - sqrt(bbb**2d0-4d0*ccc))/2d0
      write(6,*)  rydberg()*1d3* eout1
      stop
      end

      real(8) function rydberg()
      rydberg=13.6058d0
      end

      subroutine aaax()   
C -----
C Utinilty routine to calculate Im(chipm) for nmbas/=2
C set nmbas!
C fname = ChiPM*.nolfc.dat
C -----------------------------
      integer(4),parameter:: nx=1000
      integer(4):: nw,nwx,ifx=0
      real(8):: q(3),jzero,jzero2,  freq(nx), dw=0.2,rrr,iii,polinta !meV
      real(8):: rchipmm(nx),ichipmm(nx),omg,rydberg,mmom,jval,rrrx
     & ,drchi(nx),freqm(nx),czero,omgw,ddd1,ddd2,jdiff,jpp,jpm,eout
     &  ,eout1,etou2,pi
     & ,jzero3,polinta1,absq, mmnorm0,bbb,ccc,meff
      character(30) fname,aaa
      character(8) ax,bx
      integer(4)::nmbas,ix,iy
      complex(8),allocatable:: x0mean(:,:,:),x0mat(:,:),x0matinv(:,:)
     &  ,x0matinv0(:,:)
      real(8),allocatable::
     &    r0mean(:,:,:),i0mean(:,:,:),eval(:),r0mean0(:,:),eee(:),
     &    freq2(:),eee2(:)
      complex(8):: odd1,odd2,phase,img=(0d0,1d0)
      character(300) :: fff

! This is for ChiPM*.nolfc.mat ! oct2006
      pi=4*atan(1d0)
c      write(6,*)"pi=",pi
c--- Get mmom
      ifchipm=1015
      open(ifchipm, file='ChiPM.log')
      read(ifchipm,*) ax, mmom,bx, mmnorm0
      write(6,*) ax,mmom
      if(bx/='svec**2=') stop 'No svec**2 in ChiPM.log'
      write(6,*) bx,mmnorm0
      close(ifchipm)

      read(5,'(a)') fff
      read(fff,*) fname,nmbas
      nw = nx
      allocate(x0mean(nmbas,nmbas,0:nw),eee(nw),
     &    r0mean(0:nw,nmbas,nmbas) ,r0mean0(nmbas,nmbas),
     & i0mean(0:nw,nmbas,nmbas),x0mat(nmbas,nmbas)
     & ,x0matinv0(nmbas,nmbas),x0matinv(nmbas,nmbas))

      ifchipm=1015
      open(ifchipm, file=fname)
! *.nolfc.mat---> Now not 1st line skip
!      read(ifchipm,*)
      do iw=1,20
        read(ifchipm,*)q, freq(iw)
     &   , ( (r0mean(iw,ix,iy),i0mean(iw,ix,iy),ix=1,nmbas),iy=1,nmbas)
        x0mean(:,:,iw)= dcmplx(r0mean(iw,:,:),i0mean(iw,:,:))
        write(6,'("read ",3f8.4,2x,f13.7,2x,255d12.4)') 
     &   q, freq(iw) !,x0mean(iw,1:nmbas,1:nmbas)
     &   , ( (x0mean(ix,iy,iw),ix=1,nmbas),iy=1,nmbas)

!         write(6,'("ttt ",2d23.15)') sum(abs(q(:))),freq(iw)
        if(sum(abs(q(:)))==0d0.and.freq(iw)==0d0) then
           ifx=1010 
!           czero = rchipmm(iw) ; print *,q,sum(abs(q(:)))
         endif
      enddo
c---  Get czero
      if(ifx==0) then
        ifchipm0=1016
        open(ifchipm0, file='jzero')
        read(ifchipm0,*) jzero
        read(ifchipm0,*) jzero2
!     &    ((r0mean0(ix,iy),ix=1,nmbas),iy=1,nmbas)
        close(ifchipm0)
      endif


C... diagonalize
      write(6,*)' goto diagonalize --- '
      rewind(ifchipm)
C *.nolfc.mat---> Now not 1st line skip
      allocate(eval(nmbas))
      do iw=1,30
        read(ifchipm,end=1010,*)q, freq(iw)
     &   , ((r0mean(iw,ix,iy),i0mean(iw,ix,iy),ix=1,nmbas),iy=1,nmbas)
        x0mat(:,:)= dcmplx(r0mean(iw,:,:),i0mean(iw,:,:))
        x0mat(:,:)= x0mat +transpose(dconjg(x0mat))
        x0matinv= 0.5d0*x0mat
        call matcinv(nmbas,x0matinv)

ccccccccccccccccccccccccccccccccccccccc
c        write(6,*) 'test xxxxxxxxxxxx'
c        if(iw==1) then
c        odd1=x0matinv(1,2)
c        odd2=x0matinv(2,1)
c        else
c        x0matinv(1,2)= odd1
c        x0matinv(2,1)= odd2
c        endif
c        x0matinv=dreal(x0matinv)
ccccccccccccccccccccccccccccccccccccccc
c         write(6,'("x0matinv ",3f8.4,2x,f13.7,2x,255d12.4)') 
c     &   q, freq(iw), ((-x0matinv(ix,iy),ix=1,nmbas),iy=1,nmbas)
        if(iw==1) then
          jpp=  x0matinv(1,1)/4d0
          jpm= abs(-x0matinv(1,2)/4d0)
          write(6,*) 'jpp=',jpp
          write(6,*) 'jpm=',jpm
        endif
c
        if(iw==1.and.ifx/=0) then
          open(ifx, file='jzero')
          jzero = jpp + abs(jpm)
          write(6,'("x0matinv11=",2f15.5)')x0matinv(1,1)
          call diagno00(nmbas,x0matinv,eval)
          write(6,'("x0matinv eval=",255f20.15)')eval
          jzero3= dreal(x0matinv(1,1)) + abs(dreal(x0matinv(1,2)))
          jzero2=maxval(eval)
          write(ifx,*) jzero
          write(ifx,*) jzero2,jzero3 !both must be the same
          close(ifx)
        endif

c        print *,' jzero2=',jzero2
        do i=1,nmbas
          x0matinv(i,i)= x0matinv(i,i) - jzero2
        enddo
ccccccccccccccc
c          call diagno00(nmbas,x0matinv,eval)
c          write(6,'("yyy: eval=",255f20.15)')eval
ccccccccccccccc
ccccccccccccccccccccccc
c         write(6,*) 'jpp-jzero=',jpp-jzero
c         write(6,*) 'jpm=',jpm
c         write(6,'("x0matinv ",3f8.4,2x,f13.7,2x,255d12.4)') 
c     &   q, freq(iw), ((x0matinv(ix,iy),ix=1,nmbas),iy=1,nmbas)
ccccccccccccccccccc


c         write(6,'("x0inv ",3f8.4,2x,f13.7,2x,255d12.4)') 
c     &   q, freq(iw), ((dreal(x0matinv(ix,iy)),ix=1,nmbas),iy=1,nmbas)
         
c         write(6,"(3f8.4,2x,f13.7,' xxmat= ',255d12.4)") 
c     &   q, freq(iw)
c     &   , ( (x0mat(ix,iy),ix=1,nmbas),iy=1,nmbas)
        if(iw==1) x0matinv0=x0matinv
        call diagno00(nmbas,x0matinv,eval)
c        eval =  -eval *rydberg()*1d3* 4d0/mmom /4d0
        write(6,"('eval  ', 3f8.4,2x,i4,2x,f13.7,2x,2f12.6)") !,' xmat= ',255d12.4)") 
     &  q, iw, rydberg()*1d3*freq(iw), (eval(ix),ix=1,nmbas)
        eee(iw) = maxval(eval) ! &   , ( (x0mat(ix,iy),ix=1,nmbas),iy=1,nmbas)
      enddo


      do iw=1,9
        drchi(iw) = (1d0/r0mean(iw+1,1,1) - 1d0/r0mean(iw,1,1))
     &   /(freq(iw+1)-freq(iw))
        freqm(iw) = 0.5d0*(freq(iw+1)+freq(iw))
        write(6,"(' f drchi=',2f23.16)")  freqm(iw),drchi(iw)
      enddo
      omg  = 0d0
      rrr  = polinta(omg, freqm,drchi(2:8),8)
      write(6,"(2f23.16,' ! drchi=dchipm^-1/domeg \approx 1/m')") 
     & rrr, 1d0/rrr

      x0matinv0=  x0matinv0*mmnorm0/mmom**2
      meff = mmom**2/(rrr*mmnorm0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C 2x2 matrix
c      close(ifmat)
! effective  1/m = mmnorm0/mmom**2 * rrr
! mmom**2/(rrr*mmnorm0)
! Ferro case
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c-----------
 1092 format(f19.15,2x,3f19.15, f10.2,a)
 1093 format(f19.15,2x,3f19.15, f10.2,d13.5,a)
      absq=sqrt(sum(q**2))

      if(eee(1) <1d-10) then
        eout=0d0
        rrrx=-1d0
        do 
        eout = eout + 1d-6/rydberg()
        rrrx  = polinta(eout, freq, eee, 30)
        if(rrrx>0) exit
        enddo
        eouts1=eout
      else
       eout= 4d0/rrr * (-jpp+jzero-abs(jpm)) 
       write(6,1092) absq, q,  mmom,
     &  rydberg()*1d3*eout 
c        eout=0d0
c        rrrx=1d0
c        do 
c        eout = eout - 1d-6/rydberg()
c        rrrx = polinta1(eout, freq, eee, 30)
cc        write(6,*)eout,rrrx
c        if(rrrx<0) exit
c        enddo
c         eout=-99/(rydberg()*1d3)
      endif
         

c-----------
      allocate(freq2(-100:100),eee2(-100:100))
      eee2(1:100)=eee(1:100)
      freq2(1:100)=freq(1:100)
      do i=2,30
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



c--------------------------------------------
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
     &  ' ! Heis Formula AF 2atom'
      jdiff= ( (jpp - jzero)**2 - jpm**2 )* mmnorm0**2
      omgw = rydberg()*1d3* 4d0/mmom* sqrt(jdiff)
      write(6,1092)
     &  sqrt(sum(q**2)),q,  omgw,
     &  ' ! jpp-jpm formula AF 2atom'
      write(6,*)

C--------------------------------------------
      write(6,*)" --- Ferro case ---"
        write(6,1092)
     &  absq, q, rydberg()*1d3*eouts1,
     &  ' !q  m  omega(meV). pole search1'
c      write(6,1093) absq,q, mmom,
c     &  rydberg()*1d3* 4d0/rrr * (-jpp+jzero-abs(jpm))
c     &  ,(jpp-jzero) !**2-jpm**2
c     & ,' !   rrr slope Ferro'

c      write(ifmat,"('JMAT: ',3d18.10, 3x, 255d18.10)") 
c     &  q, ((x0matinv0(ix,iy),ix=1,nmbas),iy=1,nmbas) 

       ifmat=1196
       open(ifmat,file="JJJ")
       iw=1
       read(ifmat,*)q 
     &   , ( (r0mean(iw,ix,iy),i0mean(iw,ix,iy),ix=1,nmbas),iy=1,nmbas)
        x0matinv0(:,:)= dcmplx(r0mean(iw,:,:),i0mean(iw,:,:))
c      write(ifmat,"('JMAT: ',3d18.10, 3x, 255d18.10)") 
c     &  q, ((x0matinv0(ix,iy),ix=1,nmbas),iy=1,nmbas) 
        print *,' rrr readin JJJ----- meff=',meff
      bbb  = x0matinv0(1,1) + x0matinv0(2,2)
      ccc  = x0matinv0(1,1)*x0matinv0(2,2)- x0matinv0(1,2)*x0matinv0(2,1) 
      eout1 = meff*(-bbb - sqrt(bbb**2d0-4d0*ccc))/2d0

      write(6,1092)
     &  absq, q, 
     &  rydberg()*1d3* eout1, ' !  Heis Formula F 2atom'

      jdiff= (-jpp+jzero-abs(jpm))
      omgw = rydberg()*1d3* 4d0/mmom* jdiff
      write(6,1092) absq,q, omgw,
     &  ' ! jpp-jpm formula F 2atom'
      write(6,*)
      write(6,"(f14.7,' !!effective moment ')") mmom**2/(rrr*mmnorm0)
      write(6,"(f14.7,' !!static    moment ')") mmom
c
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
