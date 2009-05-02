C make interpolation ---spectrum function for metal
C   usage:   echo ChiPM*.mat |calj_interp (but do ChiPM0001.nlfc.mat, which is for q=0 first)
C   You get interpolation
C      See the sctipt calj_nlfc_metal
C
      implicit none

      integer(4),parameter:: nx=200
      integer(4):: ifx=1011,ifchipm,iw,i,istat
      real(8):: q(3),jzero,jzero2,  freq(nx), dw=0.2,rrr,iii,polinta !meV
      real(8):: rchipmm(nx),ichipmm(nx),omg,rydberg,mmom,jval,rrrx
     &,drchi(nx),freqm(nx),czero,omgw,ddd1,ddd2,jdiff,jpp,jpm,eout
     & ,eout1,etou2,pi,mmomq,sumx,cccxmx,iiix,cccx,emax,elimit
     &,jzero3,polinta1,absq, mmnorm0,bbb,ccc,meff,domega,eoutini,eoutend
      character(30) fname,aaa
      character(8) ax,bx
      integer(4)::nmbas,ix,iy
      complex(8),allocatable:: x0mean(:,:,:),x0mat(:,:),x0matinv(:,:)
     &  ,x0matinv0(:,:),xinv(:,:,:)
      real(8),allocatable::
     & r0mean(:,:,:),i0mean(:,:,:),eval(:),r0mean0(:,:),eee(:),
     & freq2(:),eee2(:),a11(:), momsite(:),usite(:),usite2(:),usite3(:),
     & uu0(:,:),uu1(:,:),uuo(:,:)
      complex(8):: odd1,odd2,phase,img=(0d0,1d0),dddc
      character(300) :: fff,ggg,hhh
      integer(4),allocatable:: imbas(:)
      real(8),allocatable:: mmnorm(:),e1(:),e2(:)
      complex(8),allocatable:: eiqrm(:)

      real(8),allocatable:: usiten(:,:),drchin(:,:)

c      read(5,*) fname
c      ifchipm=1015
c      open(ifchipm, file=fname)
c      read(ifchipm,*)
      read(5,'(a)') fff
      ggg= adjustl( fff( 1:index(fff,' '))  )
      hhh= adjustl( fff( index(fff,' '):  )  ) !this is just for older version.
      if(len(trim(hhh))/=0) then
        write(6,"(a,a)")' data file = ',trim(hhh)
        write(6,"(a,a)")' head file (for backword compatibility)) = ',trim(hhh)
      else
        write(6,"(a,a)")' data file = ',trim(ggg)
      endif
      read(ggg,*) fname
      ifchipm=1015
      open(ifchipm, file=fname)
      read(ifchipm,*) nmbas   ! number of magnetic atoms.
      allocate(imbas(nmbas),momsite(nmbas),mmnorm(nmbas),eiqrm(nmbas),
     &  e1(nmbas),e2(nmbas))
      read(ifchipm,*) imbas(1:nmbas)    ! magnetic atom index
      read(ifchipm,*) momsite(1:nmbas)  ! magnetic moment m
      read(ifchipm,*) mmnorm(1:nmbas)   ! |m|=sqrt(<m|m>) ; this is not m.
      read(ifchipm,*) (e1(ix),e2(ix),ix=1,nmbas) 
      eiqrm(1:nmbas) = dcmplx(e1(1:nmbas),e2(1:nmbas)) ! <e(iqr)|m>.  normalization?
c      write(6,*)eiqrm(1:nmbas)

      if(len(trim(hhh))/=0) then
        read(hhh,*) fname
        ifchipm=1016
        open(ifchipm, file=fname)
      endif
c----------------------------
      mmnorm0= mmnorm(1)**2
      mmom   = momsite(1)

      allocate(x0mean(nmbas,nmbas,0:nx),eee(nx),a11(0:nx),
     &    r0mean(0:nx,nmbas,nmbas) ,r0mean0(nmbas,nmbas),
     & i0mean(0:nx,nmbas,nmbas),x0mat(nmbas,nmbas)
     & ,x0matinv0(nmbas,nmbas),x0matinv(nmbas,nmbas),
     &  usite(nmbas),usite2(nmbas),usite3(nmbas)
     & ,uu0(nmbas,nmbas),uu1(nmbas,nmbas),xinv(nmbas,nmbas,nx)
     & ,uuo(nmbas,nx)  )
       uu1=0d0; uuo=0d0

       allocate(usiten(nx,nmbas),drchin(nx,nmbas))

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
ccccccccccccccccccccccccccccc

      do iw=1,nx
        read(ifchipm,*)q, freq(iw)
     &   , ( (r0mean(iw,ix,iy),i0mean(iw,ix,iy),ix=1,nmbas),iy=1,nmbas)
        x0mean(:,:,iw)= dcmplx(r0mean(iw,:,:),i0mean(iw,:,:))

c !check write
c        if(iw<4.or.mod(iw,10)==1)  then
c         write(6,'("read ",3f8.4,2x,f13.7,2x,255d12.4)') 
c     &   q, freq(iw) !,x0mean(iw,1:nmbas,1:nmbas)
c         do iy=1,nmbas
c          write(6,'($,"      ",24x,2x,255d12.4)') 
c     &    (dreal(x0mean(ix,iy,iw)),ix=1,nmbas)
c          write(6,'(" | ",255d12.4)') 
c     &    (dimag(x0mean(ix,iy,iw)) ,ix=1,nmbas)
c         enddo
c        endif

C--- x0mean is now normaized as
C  x0mean = <m|x0|m'>/sqrt(<m|m>) /sqrt(<m'|m'>)
C  where m({\bf r)} means spin density at each site.
        call matcinv(nmbas, x0mean(:,:,iw))

C ... determine usite and uu0 matrix. See kotani's paper.
        if(sum(abs(q(:)))==0d0) then
         if(iw==1) then
          ifx=1010 
          usite(1:nmbas)=
     &    matmul(dreal(x0mean(1:nmbas,1:nmbas,iw)),mmnorm(1:nmbas))
          write(6,'(" usite = ",255d12.4)')  usite(1:nmbas)

          uu0=0d0
          do i=1,nmbas
            if(abs(momsite(i))>1d-3 ) uu0(i,i)=usite(i)/mmnorm(i)
          enddo
          write(6,*)
          do iy=1,nmbas
           write(6,'("      ",24x,2x,255d12.4)') 
     &     (dreal(x0mean(ix,iy,iw)) - uu0(ix,iy),ix=1,nmbas)
          enddo
c         write(6,'(" usite= ",255d12.4)')  (usite(ix) ,ix=1,nmbas)
          write(6,'(" uu0 = ",255d12.4)') 
     &     (uu0(ix,ix) ,ix=1,nmbas)
         endif 

C ... determine uu1 matrix.
c         if(iw<=10) then
          usiten(iw,1:nmbas) = 
     &    matmul(dreal(x0mean(1:nmbas,1:nmbas,iw)),mmnorm(1:nmbas))
c          write(6,"(' usiten =',12f23.16)")  usiten(iw,:)
c          do i=1,nmbas
c          usiten(iw,i) = dreal(x0mean(i,i,iw))*mmnorm(i)
c          enddo
c         do iy=1,nmbas
c          write(6,'("vvv      ",24x,2x,255d12.4)') 
c     &    (dreal(x0mean(ix,iy,iw)),ix=1,nmbas)
c         endif
         if(iw==10) then
          do ix=1,9
            drchin(ix,:) 
     &       = (usiten(ix+1,:) - usiten(ix,:))/(freq(ix+1)-freq(ix))  !Feb24-2007
            freqm(ix) = 0.5d0*(freq(ix+1)+freq(ix))
            write(6,"(' f drchi=',12f23.16)")  freqm(ix),drchin(ix,:)
          enddo
          do i=1,nmbas
            omg  = 0d0
            rrr  = polinta(omg, freqm(2:8),drchin(2:8,i),4)
            if(abs(momsite(i))>1d-3 ) then
              uu1(i,i)= ( - momsite(i)/mmnorm(i) + rrr) /mmnorm(i)
           endif
          enddo
          write(6,'(" uu1 = ",255d12.4)') (uu1(ix,ix), ix=1,nmbas)
          write(6,'(" mom/mmnorm**2 = ",255d12.4)') 
     &     (momsite(ix)/mmnorm(ix)**2,ix=1,nmbas)
         endif 
        endif

          omg =   freq(iw)     - freq(1)
          do i=1,nmbas
            rrr =  (usiten(iw,i) - usiten(1,i)) 
            if(abs(momsite(i))>1d-3 ) then
              uuo(i,iw)= ( - omg*momsite(i)/mmnorm(i) + rrr) /mmnorm(i)
           endif
          enddo
c          write(6,'(" kkk: uuo= ",i3,255d12.4)') iw,omg,(uuo(ix,iw), ix=1,nmbas)
c          write(6,'(" kk : mom/mmnorm**2 = ",255d12.4)') 
c     &     (momsite(ix)/mmnorm(ix)**2,ix=1,nmbas)

      enddo
C .. write uu0uu1 (when q=0),  or read uu0uu1
      open(ifx,file='uu0uu1')
      if(ifx==1010) then
        write(ifx,'(255d23.15)') (uu0(ix,ix),ix=1,nmbas)
        write(ifx,'(255d23.15)') (uu1(ix,ix),ix=1,nmbas)
        do iw=1,nx
        write(ifx,'(255d23.15)') (uuo(ix,iw),ix=1,nmbas)
        enddo
      elseif(ifx==1011) then
         print *,' read uu0 uu1'
        uu0=0d0; uu1=0d0
        read(ifx,*) (uu0(ix,ix),ix=1,nmbas)
        read(ifx,*) (uu1(ix,ix),ix=1,nmbas)
        do iw=1,nx
          read(ifx,'(255d23.15)') (uuo(ix,iw),ix=1,nmbas)
        enddo
      endif
      close(ifx)
      print *,' use uu0 =',(uu0(ix,ix),ix=1,nmbas)
      print *,' use uu1 =',(uu1(ix,ix),ix=1,nmbas)

C --- get xinv = x
c      allocate(eval(nmbas))
      do iw=1,nx
c        xinv(:,:,iw) = x0mean(:,:,iw) - uu0(1:nmbas,1:nmbas) 
c     &   - freq(iw) * uu1(1:nmbas,1:nmbas) 
        xinv(:,:,iw) = x0mean(:,:,iw)  
        do ix=1,nmbas
          xinv(ix,ix,iw) = xinv(ix,ix,iw) -  uu0(ix,ix)-uuo(ix,iw) 
        enddo
c        xinv(:,:,iw) = xinv(:,:,iw) + transpose(dconjg(xinv(:,:,iw)))
c        xinv(:,:,iw) = .5d0*xinv(:,:,iw)
c        call diagno00(nmbas,xinv(:,:,iw),eval)
c        write(6,"('eval  ', 3f8.4,2x,i4,2x,f13.7,2x,255f14.8)") 
c     &  q, iw, rydberg()*1d3*freq(iw), (eval(ix),ix=1,nmbas)
c        eee(iw) = maxval(eval) ! &   , ( (x0mat(ix,iy),ix=1,nmbas),iy=1,nmbas)

c        if(iw==1) then
c          jpp=  x0mean(1,1,iw)/4d0
c          jpm= abs(-x0mean(1,2,iw)/4d0)
c          write(6,*) 'jpp=',jpp
c          write(6,*) 'jpm=',jpm
c          x0matinv0= xinv(:,:,iw)
c          jzero=1d99
c          if(sum(abs(q(:)))==0d0) jzero = jpp + abs(jpm)
c       endif
        a11(iw) = xinv(1,1,iw)
      enddo

C Get effective mass Meff--- What this means?
      do iw=1,9
        drchi(iw) = (a11(iw+1) - a11(iw))/(freq(iw+1)-freq(iw))  !Feb24-2007
        freqm(iw) = 0.5d0*(freq(iw+1)+freq(iw))
        write(6,"(' f drchi=',2f23.16)")  freqm(iw),drchi(iw)
      enddo
      omg  = 0d0
      rrr  = polinta(omg, freqm(2:8),drchi(2:8),4)
      write(6,"(2f23.16,' ! drchi=dchipm^-1/domeg \approx 1/m')") 
     & rrr, 1d0/rrr
      x0matinv0=  x0matinv0*mmnorm0/mmom**2
      meff = mmom**2/(rrr*mmnorm0)

C-------------------------------------------
      elimit=1000 !meV
      
      do iw=1,nx
        x0mat = xinv(:,:,iw)
        do i=1,nmbas
          x0mat(i,i) = x0mat(i,i) + img*1d-30 ! to avoid inversion error.
        enddo
        call matcinv(nmbas,x0mat) !this is full x0_+-  !Matrix inversion.
        dddc = sum( eiqrm*matmul(x0mat,dconjg(eiqrm)))
        rchipmm(iw) = dreal(1d0/dddc )
        ichipmm(iw) = dimag(1d0/dddc )
        write(6,"('chipm: ',3d14.6)") freq(iw),rchipmm(iw),ichipmm(iw)
      enddo

      eout= 0d0
      rrrx=-1d0
      sumx=0d0
      cccxmx=0d0
      open(106,file=trim(fname)//'.ImChiPM')
      write(106,"(' --- meV Im(Chipm) -----')")
      do 
         eout  = eout + 1d-5/rydberg()
         if(rydberg()*1d3*eout >elimit) exit
         rrrx  = polinta(eout, freq, rchipmm, nx)
         iiix  = polinta(eout, freq, ichipmm, nx)
         cccx  = -dimag(1d0/(rrrx + img*iiix))
         if(cccxmx <cccx) then
           cccxmx=cccx
           emax=eout
         endif
         if( eout> freq(nx) ) stop 'enlarge nx'
         write(106,"( f13.5, d15.7)")
     &      rydberg()*1d3*eout, cccx
         sumx=sumx+ cccx
      enddo
      close(106)
c     write(6,"('sum cccxmx=',2d13.6)") sumx,cccxmx

      istat=0
      eout=0d0
      do 
         eout  = eout + 1d-5/rydberg()
         if(rydberg()*1d3*eout >elimit) exit
         rrrx  = polinta(eout, freq, rchipmm, nx)
         iiix  = polinta(eout, freq, ichipmm, nx)
         cccx  = -dimag(1d0/(rrrx + img*iiix))
         if(istat==0.and.cccx>0.5*cccxmx) then
           eoutini=eout
           istat=1
         endif
         if(istat==1.and.cccx<0.5*cccxmx) then
            eoutend=eout
            exit
         endif
      enddo
      absq=sqrt(sum(q**2))
      write(6,*)
      write(6,"(' ---- peak of Im ChiPM ------------------')")
      write(6,"('         |q|,  q(3),  low, center, up up-low (eV),  Meff')")
      write(6,"(' FWHM: ',f19.15,2x,3f19.15,2x,3f10.3,f10.3,2x,f10.3)")
     &  absq, q, rydberg()*1d3*eoutini, rydberg()*1d3*emax, rydberg()*1d3*eoutend,
     &  rydberg()*1d3*(eoutend-eoutini), meff
      end

      real(8) function rydberg()
      rydberg=13.6058d0
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
