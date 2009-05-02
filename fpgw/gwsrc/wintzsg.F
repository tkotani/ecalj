      complex(8) function wintzsg_npm(npm,v,v0,x,wt,a,expa,we, 
     i                                nx,esmr)
C- gaussian smeared-pole version of wintz.
C takao developed from wintz by ferdi.
C Assume that each eigenvale, and w-e as a result, has the width of esmr( Ry).
c
      implicit none
      integer(4)::npm
      integer(4) :: nx,i
      real(8) :: x(nx),wt(nx),expa(nx),a,we,we2,weh,wel,weh2,wel2,cons,
     & omd,omd2,rup,rdn,sss,sig2,omd2w
      complex(8) ::v(npm*nx),v0,sum,   wintzsg1  !, wintzsg2,sum2
      complex(8) :: wwz,wintz_npm, img=(0d0,1d0)
      real(8)    :: esmr,wintsf,errsum=0d0,derfcx
      integer(4) :: ie,nav = 2000
      real(8)    :: pi=3.1415926535897932d0, rmax=2d0
      !rmax =2 is by ferdi. Is it optimum? See wintz
      real(8)::sig, smxowx, ee, omg, ww,cons1,cons2,xx,aw,eee,aw2
      integer(4) ix,verbose
      logical :: test
c      real(8)::sumgauss,sumgauss1,sumgauss2 ! I found this bug 9Aug2004
      complex(8)::sumgauss,sumgauss1,sumgauss2
      integer(4):: version

      logical ::  test_symmetric_W

c---------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' wintzsg: start!'
c      sig= 1d0 ! sig is in Ry.
c      ee = 2d0
c      do ix=1,2
c        omg= 0.1d0*(ix-1)
c        write(6,"(' omg smear=',f8.3, 2d23.15)")
c     &  omg, smxowx(omg,ee,sig), ee/(omg**2+ee**2)
c! smeared ee/(ww^2 + ee^2)
c      enddo
c      stop ' wintzsg: test end----------'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      we=1d0
c      do  i = 1,50
c        xx   = 1d-6*2**i
c        print *, ' i=',i,xx
c        sig2 =0.01d0
c        omd2w= sig2 *xx
c        cons1=(1d0 - exp (- omd2w/sig2))*we/omd2w
c        cons2 = we*( 1d0/sig2 - omd2w/sig2**2/2d0
c     &               + omd2w**2/sig2**3/6d0 - omd2w**3/sig2**4/24d0
c     &               + omd2w**4/sig2**5/120d0
c     &               - omd2w**5/sig2**6/720d0 )
c        write(6,"(' cons1=',3d25.16)")xx,cons1
c        write(6,"(' cons2=',3d25.16)")xx,cons2, cons2-cons1
c        write(6,*)
c      enddo
c      stop 'test end'


c- check part by direct integration by wintz.
c- this procedure is equivalent to whole of this routine. But time-consuming.
c      wwz =0d0
c      do ie = 1, nav
c        wex = we + 0.5d0*esmr*(-0.5d0 + dble(ie-1)/(nav-1))
c        www=1d0
c        if(ie==1.or.ie==nav) www=0.5d0
c        wwz = wwz+ www*wintz(v,v0,x,wt,a,expa,wex,nx)
cxx       write(1112,"(4d24.16)") wex,wintz(v,v0,x,wt,a,expa,wex,nx)
c      enddo
c      wintzsg1=wwz/(nav-1)
c      wintzsg=wintzsg1
c      return

      if(esmr==0d0) then
        wintzsg_npm = wintz_npm(npm,v,v0,x,wt,a,expa,we,nx)
        return
      endif
      sig  = .5d0*esmr
      sig2 = 2d0*sig**2
      we2  = we**2
      sum  = 0d0

C --- simple integration scheme.
      do  i = 1,nx
        omd   = 1d0/x(i) - 1d0
        omd2w = omd**2 + we2
! pole weight is given by   1- exp (-R^2/sig2) = 1/N \int_0^R exp(- x^2/sig2) 2 pi r dr
! Gauss theorem---> but is this correct? Not three dimentional...
        if(omd2w/sig2 > 5d-3) then
          cons = (1d0 - exp (- omd2w/sig2))/omd2w
        else
          cons = ( 1d0/sig2 - omd2w/sig2**2/2d0
     &               + omd2w**2/sig2**3/6d0  - omd2w**3/sig2**4/24d0
     &               + omd2w**4/sig2**5/120d0- omd2w**5/sig2**6/720d0 )
        endif
        sum  = sum + we*cons*( v(i)-v0*expa(i) ) *wt(i)/(x(i)**2)
        if(test_symmetric_W()) then
        else
          if(npm==2) then !Asymmetric contribution for
            sum  = sum - cons* v(i+nx)*omd*wt(i)/(x(i)**2)
          endif
        endif
      enddo

C --- Gaussian part -------------------------------------
      if(version()==-1) then !       print *,' version=',version()
        sumgauss=
     &          - .5d0*v0*dexp(we2*a**2)
     &            *( derfcx(a*we) - derfcx(sqrt(a**2+1d0/sig2)*we) )
      else
        aw = abs(a*we)
        aw2 = aw**2
        if(aw<6d0) then
          eee = we**2/sig2
          sumgauss = dsign(1d0,we)*.5d0*v0*exp(aw2)*( derfcx(sqrt(aw2+eee)) -derfcx(aw) )
c        sumgauss1=sumgauss
        else !April2004
          sumgauss =  -.5d0*v0* (2d0/(a*we)/pi)*
     &    sqrt(pi)/2d0*
     &       (1d0 - 1d0/2d0/aw2*
     &         (1d0 - 3d0/2d0/aw2*
     &           (1d0 - 5d0/2d0/aw2*
     &             (1d0 - 7d0/2d0/aw2* 
     &               (1d0 - 9d0/2d0/aw2*
     &                 (1d0 - 11d0/2d0/aw2 
     &     ))))))
c       sumgauss2=sumgauss
        endif
      endif
      wintzsg_npm = -sum/pi  + sumgauss
      if(verbose()>90) 
     & write(1116, "(' we sig wintzsg_npm =',2f8.3,4f14.6)" )we,sig,wintzsg_npm
      return
c---test code
      print *,' sum= ', sum
      print *,' awe1= ', a*we
      print *,' erf1= ', derfcx(a*we)
      print *,' awe2= ', sqrt(a**2+1d0/sig2)*we
      print *,' erf2= ', derfcx(sqrt(a**2+1d0/sig2)*we)
      print *,' dexp= ', dexp(we2*a**2)
      print *,' erf3= ', derfcx(aw)
      print *,' erf4= ', derfcx(sqrt(aw2+eee))
      print *,' aw2= ', aw2
      aw = a*we
      write(6, "(' wintzsg_npm chk1=',2f8.3, 4d18.10)" )we,sig,
     &  - .5d0*v0*  dexp(aw**2)*( derfcx(aw) - derfcx(sqrt(a**2+1d0/sig2)*we)  )
      write(6, "(' wintzsg_npm chk2=',2f8.3, 4d18.10)" )we,sig,sumgauss1
      write(6, "(' wintzsg_npm chk3=',2f8.3, 4d18.10)" )we,sig,sumgauss2
      return
      end




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--------------------------------------------------------------------------
      real(8) function smxowx(ww,ee,sig)
      implicit none
      real(8),parameter:: pi=3.1415926535897932d0
      real(8):: sq2pi,s2int
      real(8):: sum, fac,wwsig,s2,ww,ee,sig,dn2fac
      integer(4)::  mxn,n
      sq2pi=sqrt(2d0*pi)
      fac= exp( - ee**2/sig**2/2d0)*ee/sig
      mxn= 30 !???
      wwsig= ww/sig
      s2  = s2int(wwsig)
      print *,'wwsig=',wwsig
      print *,' s2=',s2
      sum = fac*s2
      dn2fac= 1d0
      do n=2,mxn
        fac = fac*(ee/sig)**2 /(2*n-1)/(2*n-2)
        dn2fac = dn2fac*(2*n-1)
        s2  = sq2pi*dn2fac - wwsig**2*s2
        sum = sum + fac*s2
        print *,' n sum',n,sum/sqrt(2*pi)
      enddo
      smxowx=sum/sqrt(2d0*pi)/sig
      end
      real(8) function s2int(beta)
      implicit none
      real(8),parameter:: pi=3.1415926535897932d0
      real(8):: derfc,sqpi,beta,derf
      real(8):: mu
      sqpi= sqrt(pi)
      mu  = 1d0/sqrt(2d0)
c      print *, '----------------beta*mu=',beta*mu
      s2int = sqpi/mu - pi*beta*exp(beta**2/2d0)*derfc(beta*mu)
      end


c -------------------------------------------------------------------
c --------------------------------------------------------------------
      complex(8) function wintzav(v,v0,x,wt,a,expa,we, 
     i                                nx,esmr)
C- Smeared-pole version of wintz. takao developed from wintz by ferdi.
C Assume thaat each eigenvale, and w-e as a result, has the width of esmr( Ry).
c
      implicit none
      integer(4) :: nx,i
      real(8) :: x(nx),wt(nx),expa(nx),a,we,we2,weh,wel,weh2,wel2,cons,
     & omd,omd2,rup,rdn,sss
      complex(8) ::v(nx),v0,sum,   wintzav1  !, wintzav2,sum2
      complex(8) :: wwz,wintz
      real(8)    :: esmr,wintsf,errsum=0d0
      integer(4) :: ie,nav = 2000
      real(8)    :: pi=3.1415926535897932d0, rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
     &
c---------------------------------------------------
c      print *,' wintzav:'
      if(esmr==0d0) then
        wintzav=wintz(v,v0,x,wt,a,expa,we,nx)
        return
      endif
c-----------------------------------------------------
c- check part by direct integration by wintz.
c- this procedure is equivalent to whole of this routine. But time-consuming.
c      wwz =0d0
c      do ie = 1, nav
c        wex = we + 0.5d0*esmr*(-0.5d0 + dble(ie-1)/(nav-1))
c        www=1d0
c        if(ie==1.or.ie==nav) www=0.5d0
c        wwz = wwz+ www*wintz(v,v0,x,wt,a,expa,wex,nx)
cxx       write(1112,"(4d24.16)") wex,wintz(v,v0,x,wt,a,expa,wex,nx)
c      enddo
c      wintzav1=wwz/(nav-1)
c      wintzav=wintzav1
c      return

c frequency integral
      we2  = we**2
      weh  = we + 0.25d0*esmr !because we is in a.u and esmr is in Ry.
      wel  = we - 0.25d0*esmr
      weh2 = weh**2
      wel2 = wel**2
      sum  = 0d0
c
c      sum2 = 0d0
c      wintzav2 = 0d0

ccase1 --------------------
      if ( max(dabs(weh),dabs(wel)) < rmax/a) then
ccase2 --------------------
ccase2 if (.false.) then
        do  i = 1,nx
          omd  = 1d0/x(i) - 1d0
          omd2 = omd**2
c   x = 1/(1+w')
c   We smear each pole with a pole distribution between wel ans weh,
c   So  it means \int_wel^weh dwe [ wt(i)/(x(i)**2) * we/(we2 + omd**2)]
          cons  = 0.5d0*log( (weh2+omd2)/(wel2+ omd2)) /(x(i)**2)
          sum   = sum + wt(i)*(v(i) - v0*expa(i))*cons
ccccccccccccccccccccccccccccccccccccccccccccccccc
c        sum2  = sum2 + wt(i)*v(i) *cons
c      write(6,"(d15.6,2x,2d15.6,2x,2d15.6)")
c     & x(i), v(i)*cons,v0*expa(i)*cons
c      write(6,"(d15.6,2x,2d16.6,2x,2d16.6)")
c     & x(i), v(i),v0*expa(i)
ccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo

        if(abs(weh)>abs(wel)) then
          sss = 1d0
          rup = abs(weh)
          rdn = abs(wel) !;    iix=1
        else
          sss = -1d0
          rup = abs(wel)
          rdn = abs(weh) !;   iix=11
        endif
        wintzav  = - sum/pi - 0.5d0*v0* sss* wintsf(a*rdn,a*rup)/a
c       wintzav2 = - sum2/pi
ccccccccccccccccccccccccccccccccccccccccccccccccc
c      xxx=0d0
c      write(6,"(d15.6,2x,2d16.6,2x,2d16.6)")
c     & xxx, v0 ,v0
c      write(6,"(' - sum/pi =',2d13.6)") -sum/pi
c      write(6,"(' wintzav  =',2d13.6)")   wintzav
c      write(6,"(' - sum2/pi =',2d13.6)") -sum2/pi
c      stop ' ==== test end of wintz ========'
ccccccccccccccccccccccccccccccccccccccccccccccccc
      else
c simple integration scheme.
        do  i = 1,nx
          omd  = 1d0/x(i) - 1d0
          omd2 = omd**2
          cons =  0.5d0*log( (weh2+omd2)/(wel2+ omd2)) /(x(i)**2)
          sum  = sum + wt(i)*v(i)*cons
        enddo
        wintzav = -sum/pi  !;   iix =2
      endif
      wintzav = wintzav/(0.5d0*esmr) !because esmr is in Ry.

cccccccccccccccccccccccccccc
c      write(1116, "(' wel weh wintzav =',2f8.3,4f14.6)" )wel,weh,wintzav
ccccccccccccccccccccccccccc

c      wintzav2 = wintzav2/(0.5d0*esmr)
c      write(1116, "(' wintzav diff =',4d13.6,f10.6)" ) wintzav,wintzav2
c     &   ,abs(wintzav2-wintzav)/abs(wintzav)
cccccccc
c      stop ' ==== test end of wintz ========'
c      if(wintzav/=0d0 ) then
c        err = abs(wintzav-wintzav1)
c        if(err>1d-6) then
c          write(3302,"(3d24.16, i3,' err= ',2d13.5)")
c     &    we, wintzav, iix, err
c        endif
c     endif
cccccccccccccccccccccccccccc
      end

c----------------------------------------------------
      real(8) function wintsf(rdn,rup)
c integral \int^rup_rdn exp(x**2)*derfc(x)
c up to x <rmax<15 which is limited by the range of derfc(x).
      implicit real*8 (a-h,o-z)
      integer(4),parameter:: nd=2001
      real(8) :: polinta, rmax=14.9d0
      real(8),save :: eeint(nd),xval(nd)
      logical :: init=.true.
c--------------
      if(rup>rmax) stop 'wintsf: too large rup'
      if(rup<rdn ) stop 'wintsf: rup<rdn'
      if(init) then
        eeint(1) = 0d0
        xval(1)  = 0d0
        do ix= 2,nd
          x0= rmax*dble(ix-2)/(nd-1)
          x = rmax*dble(ix-1)/(nd-1)
          xm= .5d0*(x + x0)
          xval(ix) = x
          sum  = 1/6d0 * exp(x0**2)*derfc(x0)
     &         + 4/6d0 * exp(xm**2)*derfc(xm)
     &         + 1/6d0 * exp(x **2)*derfc(x )
          eeint(ix) = eeint(ix-1) + rmax/(nd-1) * sum
ccccccccccccccccccccccc
c        write(1113,"(4d24.16)") x,derfc(x),exp(x**2)*derfc(x),eeint(ix)
ccccccccccccccccccccccc
        enddo
        init=.false.
      endif
      wintsf = polinta(rup,xval,eeint,nd) - polinta(rdn,xval,eeint,nd)
      end


c--------------------------------------------------------------------
c--------------------------------------------------------------------
      complex(8) function wintz_npm(npm,v,v0,x,wt,a,expa,we,
     d                                nx)
c takao complex version of wint by ferdi
c 92.02.24
c wintz = < [w'=-inf,inf] (i/2pi) v(w')/(w+w'-e) >

c the integration is done along the imaginary axis
c w' ==> iw', w' is now real
c wintz = - < [w'=0,inf] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2} >

c transform: x = 1/(1+w')
c wintz = - < [x=0,1] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >

c the integrand is peak around w'=0 or x=1 when we=w-e=0
c to handel the problem, add and substract the singular part as follows:
c wintz = - < [x=0,1] {v(iw') - v(0)exp(-a^2 w'^2)}
c                   *(1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
c        - (1/2) v(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
c the second term of the integral can be done analytically, which
c results in the last term

c when we=w-e ==> 0, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w')
c the integral becomes -v(0)/2

c v       = v(iw')
c v0      = v(0)
c x       = s.o.
c wt      = weights for integration
c a       = a constant determining the range of exp(-a^2 w'^2)
c expa(x) = exp(-a^2 w'^2)
c we      = w - e
c nx      = no. x points

      implicit real*8 (a-h,o-z)
      dimension x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx),v0,sum,img=(0d0,1d0),wintz
      real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
      data pi/3.1415926535897932d0/,
cmod 93.06.23
c    d     tol/1.d-6/
     d     tol/1.d-8/

c if w = e the integral = -v(0)/2
cmod 93.06.23
c     if (dabs(we) .lt. tol) then
c     wintz       = -0.5d0*v0
c     return
c     endif
      if (dabs(we) .lt. tol) then
        print *, ' we=',we
        stop 'wintz: |w-e| < tol'
      endif
c frequency integral
      we2        = we*we
      sum        = 0.d0

      if (dabs(we) < rmax/a) then
        do       i = 1,nx
          omd   = 1d0/x(i) - 1d0
          onemx      = 1.d0 - x(i)
          cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + we*cons*(v(i) - v0*expa(i))*wt(i)
          if(npm==2) then !Asymmetric contribution for
            sum  = sum - cons* v(i+nx)*omd*wt(i)
          endif
        enddo

        wintz      = -sum/pi 
     .            - 0.5d0*v0*dsign(1.d0,we)*dexp(we2*a*a)
     .                                     *derfc(a*dabs(we))
      else
        do       i = 1,nx
          omd   = 1d0/x(i) - 1d0 !this was missing. I added this at 25June2008
          onemx      = 1.d0 - x(i)
          cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + we*cons*v(i)*wt(i)
          if(npm==2) then !Asymmetric contribution for
            sum  = sum - cons* v(i+nx)*omd*wt(i)
          endif
        enddo

        wintz      = -sum/pi
      endif

      wintz_npm= wintz
      return
      end

c-------------------------------------------------------------
      complex(8) function wintz(v,v0,x,wt,a,expa,we,
     d                                nx)
c takao complex version of wint by ferdi
c 92.02.24
c wintz = < [w'=-inf,inf] (i/2pi) v(w')/(w+w'-e) >

c the integration is done along the imaginary axis
c w' ==> iw', w' is now real
c wintz = - < [w'=0,inf] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2} >

c transform: x = 1/(1+w')
c wintz = - < [x=0,1] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >

c the integrand is peak around w'=0 or x=1 when we=w-e=0
c to handel the problem, add and substract the singular part as follows:
c wintz = - < [x=0,1] {v(iw') - v(0)exp(-a^2 w'^2)}
c                   *(1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
c        - (1/2) v(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
c the second term of the integral can be done analytically, which
c results in the last term

c when we=w-e ==> 0, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w')
c the integral becomes -v(0)/2

c v       = v(iw')
c v0      = v(0)
c x       = s.o.
c wt      = weights for integration
c a       = a constant determining the range of exp(-a^2 w'^2)
c expa(x) = exp(-a^2 w'^2)
c we      = w - e
c nx      = no. x points

      implicit real*8 (a-h,o-z)
      dimension x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx),v0,sum
      real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
      data pi/3.1415926535897932d0/,
cmod 93.06.23
c    d     tol/1.d-6/
     d     tol/1.d-8/

c if w = e the integral = -v(0)/2
cmod 93.06.23
c     if (dabs(we) .lt. tol) then
c     wintz       = -0.5d0*v0
c     return
c     endif
      if (dabs(we) .lt. tol) then
        print *, ' we=',we
        stop 'wintz: |w-e| < tol'
      endif
c frequency integral
      we2        = we*we
      sum        = 0.d0

      if (dabs(we) < rmax/a) then
        do       i = 1,nx
          onemx      = 1.d0 - x(i)
          cons       = we/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + cons*(v(i) - v0*expa(i))*wt(i)
        enddo

        wintz      = -sum/pi 
     .           - 0.5d0*v0*dsign(1.d0,we)*dexp(we2*a*a)
     .                                    *derfc(a*dabs(we))
      else
        do       i = 1,nx
          onemx      = 1.d0 - x(i)
          cons       = we/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + cons*v(i)*wt(i)
        enddo

        wintz      = -sum/pi
      endif

      return
      end

c--------------------------------------------------------------------
      complex(8) function alagr3z_even (x,xi,fi)
      implicit real*8 (a-h,o-z)
c evne function version of alagr3z
c return the interpolated value on x for fi(xi).
c We assume even function.
c Imag part is corrected to be >0
      dimension  xi(3)
      complex(8) :: fi(3) ,alagr3z
      real(8) :: amatinv(3,3),amat(3,3),ratio

c--- Assure that dimag(fi) is positive.
      if(x<0d0) stop ' alagr3z: x<0d0'
      do i=1,3
        if(dimag(fi(i))>1d-8) then
          ratio= abs(dimag(fi(i))/dreal(fi(i)))
          if( ratio >1d-6) then !jan2004. I added a bit of alloance
            print *,' alagr3z: error! dimag(fi(i)) is positive',i,xi,fi
            stop    ' alagr3z: error! dimag(fi(i)) is positive'
          endif
        endif
      enddo
c
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
      call dinv33(amat,0,amatinv,detxx)
      alagr3z = dcmplx (
     &          sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ),
     &          sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )

c--- imag part of alagr3z is positive because interpolation could give
c    negative values even if all fi(i) <0.
c
      if(dimag(alagr3z)>0d0 ) then
        alagr3z = dcmplx( dreal(alagr3z),0d0)
      endif
      alagr3z_even=alagr3z
      end

c
c-------------------------------------
      subroutine matzwz(zw,zmel, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,ngb
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
      complex(8), allocatable :: CC(:,:,:)
      allocate(CC(ngb,nstate,ntp0) )
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntp0)
      do itp = 1,ntp0
        do  it = 1,nstate
          zwz(it,itp) = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itp))
        enddo
      enddo
      deallocate(CC)
      end
c-------------------------------
      subroutine gen_ua(ratio,niw,freqx, expa_,ua_)
      implicit none
      integer(4):: ix,niw,verbose
      real(8):: ratio,freqw, freqx(niw), expa_(niw),ua_  !, ua_min=1d0
      ix=niw
      freqw = (1d0 - freqx(ix))/ freqx(ix)

c      if( log(ratio) > -1d0/freqw/ua_min**2) then
c         print *,' ratio =',ratio
c         print 'gen_ua: ratio >1 strange ratio...'
c         stop  'gen_ua: ratio >1 strange ratio...'
c        ua_ = ua_min
      if( ratio > 1d0 ) then
        print *, ' ratio =',ratio
        print *, ' gen_ua: ratio >1 strange ratio...'
        stop  ' gen_ua: ratio >1 strange ratio...'
      else
        ua_ = sqrt(-1d0/freqw*log(ratio))
      endif
      do ix = 1,niw
        freqw     = (1d0 - freqx(ix))/ freqx(ix)
        expa_(ix) = exp(-(ua_*freqw)**2)
      enddo
      end
      subroutine gen_uaa(ratio,freqx, ua_)
      implicit none
      integer(4):: ix,niw,verbose
      real(8):: ratio,freqw, freqx, ua_
      freqw = (1d0 - freqx)/ freqx
      if( ratio > 1d0 ) then
        print *, ' ratio =',ratio
        print *, ' gen_uaa: ratio >1 strange ratio...'
        stop  ' gen_uaa: ratio >1 strange ratio...'
      else
        ua_ = sqrt(-1d0/freqw*log(ratio))
      endif
      end

      subroutine gen_expa(niw,freqx,ua_,  expa_)
      implicit none
      integer(4):: ix,niw
      real(8)::  freqw, freqx(niw), expa_(niw), ua_
      do ix = 1,niw
        freqw     = (1d0 - freqx(ix))/ freqx(ix)
        expa_(ix) = exp(-(ua_*freqw)**2)
      enddo
      end

c-------------------------------------------------------------------
      subroutine matzwz2(iSigma_en ,zw,zmel, ntq,nstate,ngb, zwz) !sf 23May02
      implicit none
      integer(4) :: nstate,ntq,itp,itpp,it,ngb,iSigma_en
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntq)
      complex(8) :: zwz(nstate,ntq,ntq) ,zdotc                   !sf 23May02
      complex(8), allocatable :: CC(:,:,:)
      allocate(CC(ngb,nstate,ntq) )
      !write(*,'("10.1 matzwz begin ",$)')
      !call cputid(0)
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntq) !most time consuming part
      !write(*,'("10.2 matzwz after matm ",$)')
      !call cputid(0)
      if (iSigma_en==1.or.iSigma_en==5) then                                       !sf 23May02
        zwz=0d0
        do itp = 1,ntq
          do  it = 1,nstate
c        zwz(it,itp,itp)=sum(dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itp))
            zwz(it,itp,itp)=
     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itp),1 )
          enddo
        enddo
      elseif (iSigma_en==2 .or. iSigma_en==3 ) then
        do itpp= 1,ntq
          do itp = 1,ntq
            do  it = 1,nstate
              zwz(it,itp,itpp)=
     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itpp),1 )
c     &    sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itpp) )      !sf 23May02
            enddo
          enddo
        enddo
      else
        stop "sxcf, matzwz2: iSigma_en /= 0,1,2,3"
      endif
      !write(*,'("10.3 matzwz end ",$)')
      !call cputid(0)

      deallocate(CC)
      end

c--------------------------------------------------------------------
      complex(8) function alagr3z2(x,xi,fi,i,j)  !sf 24 May02
      implicit real*8 (a-h,o-z)
c evne function version of alagr3z
c return the interpolated value on x for fi(xi).
c We assume even function.
c Imag part is corrected to be >0
      integer(4) :: i,j,n
      dimension  xi(3)
      complex(8) :: fi(3)
      real(8) :: amatinv(3,3),amat(3,3),ratio

c--- Assure that dimag(fi) is positive.
      if(x<0d0) stop ' alagr3z2: x<0d0'
      if (i==j) then  !sf 24 May02 only diagonal Im[Sigma_ii] should be negative
        do n=1,3
          if(dimag(fi(n))>1d-8) then
            ratio= abs(dimag(fi(i))/dreal(fi(i)))
            if( ratio >1d-8) then !jan2004. I added a bit of alloance
              print *,' alagr3z: error! dimag(fi(i)) is positive',i,xi,fi
              stop    ' alagr3z: error! dimag(fi(i)) is positive'
            endif
          endif
        enddo
      endif !(i=j)
c
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
      call dinv33(amat,0,amatinv,detxx)
      alagr3z2 = dcmplx (
     &          sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ),
     &          sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )

c--- imag part of alagr3z is positive because interpolation could give
c    negative values even if all fi(i) <0.
c
      if(i==j) then
        if(dimag(alagr3z2)>0d0 ) then
          alagr3z2 = dcmplx( dreal(alagr3z2),0d0)
        endif
      endif !(i==j)
      end

      subroutine timeshow(info)
      character*(*) :: info
      write(*,'(a,$)')info
      call cputid(0)
      end




c======================================================================
      subroutine matzwz_notusedx(zw,zmel, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,ngb
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
c     & ,wx(ngb)
c --- This routine is really time consuming! ---
      do itp = 1,ntp0
        do  it = 1,nstate
          zwz(it,itp) = sum( 
     &   dconjg(zmel(:,it,itp))*matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)))
c           wx = matmul(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp)) * wx )
        enddo
      enddo
      end
c---------------------------------------------
      subroutine matzwzs_notused(zw,zmel, ntq,nstate,ngb, zwz) !sf 23May02
c I have not checked it yet...
      implicit none
      integer(4) :: nstate,ntq,itp,itpp,it,ngb,iSigma_en,igb,igb2
      complex(8) :: zw(ngb,ngb),zwx(ngb,ngb),zmel(ngb,nstate,ntq)
      complex(8) :: zwz(nstate,ntq,ntq) ,zdotc,zz3(ntq),zz2
      complex(8), allocatable :: CC(:,:,:),zwzx(:,:,:),zmel1_(:,:,:)
      allocate(zwzx(nstate,ntq,ntq),zmel1_(ntq,ngb,nstate))
      do itpp= 1,ntq
        do it  = 1,nstate
          zmel1_(itpp,1:ngb,it) = zmel(1:ngb,it,itpp)
        enddo
      enddo
      zwx = zw
      do igb=1,ngb
        zwx(igb,igb)=.5d0*zwx(igb,igb)
      enddo
      do itp=1,ntq
        do it =1,nstate
c        zmel1(:)=dconjg(zmel(:,it,itp))
          zz3=0d0
          do igb2=1,ngb
            zz2 = zdotc(igb2, zmel(1,it,itp),1,zwx(1,igb2),1 )
            call zaxpy(ntq,zz2,zmel1_(1,igb2,it),1,zz3,1)
          enddo
          zwzx(it,itp,:)= zz3(:)
        enddo
      enddo
      do itpp=1,ntq
        do itp =1,ntq
          do it  =1,nstate
            zwz(it,itp,itpp)=zwzx(it,itp,itpp)+dconjg(zwzx(it,itpp,itp))
          enddo
        enddo
      enddo
      deallocate(zwzx,zmel1_)
      return
c-----------------------------------------------
      Allocate(CC(ngb,nstate,ntq) )
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntq) !most time consuming part
      !write(*,'("10.2 matzwz after matm ",$)')
      !call cputid(0)
      do itpp= 1,ntq
        do itp = 1,ntq
          do  it = 1,nstate
            zwz(it,itp,itpp)=
     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itpp),1 )
          enddo
        enddo
      enddo
      deallocate(CC)
      end
c      subroutine matzwz(zw,zmel, ntq,nstate,ngb, zwz)
c      implicit none
c      integer(4) :: nstate,ntq,itp,it,ngb
c      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntq),zwz(nstate,ntq)
c      complex(8), allocatable :: CC(:,:,:)
c      allocate(CC(ngb,nstate,ntq) )
c
c      call matm(zw,zmel,cc, ngb, ngb, nstate*ntq)
c
c      do itp = 1,ntq
c      do  it = 1,nstate
c        zwz(it,itp) = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itp))
c      enddo
c      enddo
c
c      deallocate(CC)
c      end

c-------------------------------------
c   This does not work.
c      subroutine matzwz(zw,zmel, ntp0,nstate,ngb, zwz)
c      implicit none
c      integer(4) :: nstate,ntp0,itp,it,ngb
c      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
c     & ,wx(ngb)
c --- This routine is really time consuming! ---
c        do itp = 1,ntp0
c        do  it = 1,nstate
c        call mml(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp), wx(1:ngb))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp))*wx )
c           wx = matmul(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp)) * wx )
c        enddo
c        enddo
c      end

c -----------------------------------
      complex(8) function alagr3z (x,xi,fi)
c complex version of fi
c 92.03.02
c three-point interpolation with arbitrary mesh
c given a function fi at three points xi, the routine interpolates
c the function at x
c f(x) = [ { (x-x2)(x-x3) } / { (x1-x2)(x1-x3) } ] f1
c      + [ { (x-x1)(x-x3) } / { (x2-x1)(x2-x3) } ] f2
c      + [ { (x-x1)(x-x2) } / { (x3-x1)(x3-x2) } ] f3
c
c x  = the point at which the function is to be interpolated
c xi(3) = points where the function is given
c fi(3) = the function at xi
      implicit real*8 (a-h,o-z)
      dimension xi(3)
      complex(8) :: fi(3)
      xx1        = x-xi(1)
      xx2        = x-xi(2)
      xx3        = x-xi(3)
      x12        = xi(1)-xi(2)
      x13        = xi(1)-xi(3)
      x23        = xi(2)-xi(3)
      alagr3z    = xx2*xx3*fi(1)/(x12*x13)
     .           - xx1*xx3*fi(2)/(x12*x23)
     .           + xx1*xx2*fi(3)/(x13*x23)
      return
      end

c--------------------------------------------------------------------
      complex(8) function alagr3zz(x,xi,fi)  !jan 2008 takao. starting from alagr3z2
c evne function version of alagr3z !
c return the interpolated value on x for fi(xi).
c We assume even function.
c Imag part is corrected to be >0
      real(8)::  xi(3), amatinv(3,3),amat(3,3),x
      complex(8) :: fi(3)
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
      call dinv33(amat,0,amatinv,detxx)
      alagr3zz = dcmplx (
     &          sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ),
     &          sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )
c      write(6,"(' coeff amatinv = ',9d11.3)") matmul(amatinv,dimag(fi))
c      write(6,"(' coeff    x    = ',9d11.3)") x**2,x**4
c--- imag part of alagr3z is positive because interpolation could give
c    negative values even if all fi(i) <0.
      if(dimag(alagr3zz)>0d0) alagr3zz = dcmplx( dreal(alagr3zz),0d0)
      end
