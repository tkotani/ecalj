      subroutine wintzsg_npm_wgtim(npm,x,wt,a,expa,we,nx,esmr, wgtim)
      implicit none
      integer(4),intent(in)::npm,nx
c      complex(8),intent(in) ::v(npm*nx),v0
      real(8),intent(in):: x(nx),wt(nx),expa(nx),we,esmr,a
      integer(4) :: i
      real(8):: we2,weh,wel,weh2,wel2,cons,
     & omd,omd2,rup,rdn,sss,sig2,omd2w
c      complex(8) :: wwz,wintz_npm, img=(0d0,1d0),sum,   wintzsg1 
      real(8)    :: wintsf,errsum=0d0,derfcx,derfc
      integer(4) :: ie,nav = 2000
      real(8),parameter :: pi=3.1415926535897932d0, rmax=2d0
      !rmax =2 is by ferdi. Is it optimum? See wintz
      real(8)::sig, smxowx, ee, omg, ww,cons1,cons2,xx,aw,eee,aw2
      integer(4) ix,verbose
      logical :: test
      real(8):: wgtim(0:npm*nx)
      if(esmr==0d0) then
        call wintz_npm_wgtim(npm,x,wt,a,expa,we,nx, wgtim)
        return
      endif
      sig  = .5d0*esmr
      sig2 = 2d0*sig**2
      we2  = we**2
      wgtim=0d0
!! simple integration scheme.
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
       wgtim(i)= wgtim(i)+ we*cons*wt(i)/(x(i)**2)*(-1d0/pi)
       wgtim(0)= wgtim(0)- we*cons*expa(i)*wt(i)/(x(i)**2)*(-1d0/pi)
c        sum  = sum + we*cons*( v(i)-v0*expa(i) ) *wt(i)/(x(i)**2)
          if(npm==2) then !Asymmetric contribution for
            wgtim(i+nx) = wgtim(i+nx) - cons*omd*wt(i)/(x(i)**2)*(-1d0/pi)
          endif
      enddo
!! --- Gaussian part -------------------------------------
        aw = abs(a*we)
        aw2 = aw**2
        if(aw<6d0) then
          eee = we**2/sig2
          wgtim(0)=wgtim(0)+ dsign(1d0,we)*.5d0*exp(aw2)*( derfcx(sqrt(aw2+eee)) -derfcx(aw) )
        else !April2004
          wgtim(0)= wgtim(0) -.5d0* (2d0/(a*we)/pi)*
     &    sqrt(pi)/2d0*
     &       (1d0 - 1d0/2d0/aw2*
     &         (1d0 - 3d0/2d0/aw2*
     &           (1d0 - 5d0/2d0/aw2*
     &             (1d0 - 7d0/2d0/aw2* 
     &               (1d0 - 9d0/2d0/aw2*
     &                 (1d0 - 11d0/2d0/aw2 
     &     ))))))
      endif
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wintz_npm_wgtim(npm,x,wt,a,expa,we,nx, wgtim)
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
      integer:: nx,i,npm
      real(8):: x(nx),wt(nx),expa(nx)
c      complex(8) ::v(nx),v0,sum,img=(0d0,1d0),wintz
      real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
      real(8),parameter:: pi=3.1415926535897932d0,tol=1d-8
      real(8):: wgtim(0:npm*nx)
c if w = e the integral = -v(0)/2
c frequency integral
      if (dabs(we)< tol) call rx1( 'wintz: |w-e| < tol',we)
      we2  = we*we
      sum  = 0d0
      wgtim = 0d0
      if (dabs(we) < rmax/a) then
         do       i = 1,nx
            omd   = 1d0/x(i) - 1d0
            onemx      = 1.d0 - x(i)
            cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
            wgtim(i)= wgtim(i)+ we*cons*wt(i)*(-1d0/pi)
            wgtim(0)= wgtim(0)+ we*cons*(-expa(i))*wt(i)*(-1d0/pi)
            if(npm==2) then     !Asymmetric contribution for
               wgtim(i+nx) = wgtim(i+nx) - cons*omd*wt(i)*(-1d0/pi)
            endif
         enddo
         wgtim(0)= wgtim(0)-0.5d0*dsign(1.d0,we)*dexp(we2*a*a)*derfc(a*dabs(we))
      else
         do       i = 1,nx
            omd   = 1d0/x(i) - 1d0 !this was missing. I added this at 25June2008
            onemx      = 1.d0 - x(i)
            cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
            wgtim(i)= wgtim(i)+ we*cons*wt(i)*(-1d0/pi)
            if(npm==2) then     !Asymmetric contribution for
               wgtim(i+nx)= wgtim(i+nx) - cons* omd*wt(i)*(-1d0/pi)
            endif
         enddo
      endif
      end
!> gaussian smeared-pole version of wintz.
!! takao developed from wintz by ferdi.
!! Assume that each eigenvale, and w-e as a result, has the width of esmr( Ry).
      complex(8) function wintzsg_npm(npm,v,v0,x,wt,a,expa,we, 
     i                                nx,esmr)
      implicit none
      integer(4),intent(in)::npm,nx
      complex(8),intent(in) ::v(npm*nx),v0
      real(8),intent(in):: x(nx),wt(nx),expa(nx),we,esmr

      integer(4) :: i
      real(8):: a,we2,weh,wel,weh2,wel2,cons,
     & omd,omd2,rup,rdn,sss,sig2,omd2w
      complex(8) :: wwz,wintz_npm, img=(0d0,1d0),sum,   wintzsg1 
      real(8)    :: wintsf,errsum=0d0,derfcx,derfc
      integer(4) :: ie,nav = 2000
      real(8)    :: pi=3.1415926535897932d0, rmax=2d0
      !rmax =2 is by ferdi. Is it optimum? See wintz
      real(8)::sig, smxowx, ee, omg, ww,cons1,cons2,xx,aw,eee,aw2
      integer(4) ix,verbose
      logical :: test
      complex(8)::sumgauss,sumgauss1,sumgauss2
      if(esmr==0d0) then
        wintzsg_npm = wintz_npm(npm,v,v0,x,wt,a,expa,we,nx)
        return
      endif
      sig  = .5d0*esmr
      sig2 = 2d0*sig**2
      we2  = we**2
      sum  = 0d0

!! simple integration scheme.
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
c        if(test_symmetric_W()) then
c        else
          if(npm==2) then !Asymmetric contribution for
            sum  = sum - cons* v(i+nx)*omd*wt(i)/(x(i)**2)
          endif
c        endif
      enddo
!! --- Gaussian part -------------------------------------
c      if(version()==-1) then !       print *,' version=',version()
c        sumgauss=
c     &          - .5d0*v0*dexp(we2*a**2)
c     &            *( derfcx(a*we) - derfcx(sqrt(a**2+1d0/sig2)*we) )
c      else
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
      s2int = sqpi/mu - pi*beta*exp(beta**2/2d0)*derfc(beta*mu)
      end
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
      if(esmr==0d0) then
        wintzav=wintz(v,v0,x,wt,a,expa,we,nx)
        return
      endif
c frequency integral
      we2  = we**2
      weh  = we + 0.25d0*esmr !because we is in a.u and esmr is in Ry.
      wel  = we - 0.25d0*esmr
      weh2 = weh**2
      wel2 = wel**2
      sum  = 0d0
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
      end

c----------------------------------------------------
      real(8) function wintsf(rdn,rup)
c integral \int^rup_rdn exp(x**2)*derfc(x)
c up to x <rmax<15 which is limited by the range of derfc(x).
      implicit none
      integer:: ix
      integer,parameter:: nd=2001
      real(8) :: polinta, rmax=14.9d0,rdn,rup,sumx,x,x0,xm,derfc
      real(8),save :: eeint(nd),xval(nd)
      logical :: init=.true.
      if(rup>rmax) call rx( 'wintsf: too large rup')
      if(rup<rdn ) call rx( 'wintsf: rup<rdn')
      if(init) then
        eeint(1) = 0d0
        xval(1)  = 0d0
        do ix= 2,nd
          x0= rmax*dble(ix-2)/(nd-1)
          x = rmax*dble(ix-1)/(nd-1)
          xm= .5d0*(x + x0)
          xval(ix) = x
          sumx  = 1/6d0 * exp(x0**2)*derfc(x0)
     &         + 4/6d0 * exp(xm**2)*derfc(xm)
     &         + 1/6d0 * exp(x **2)*derfc(x )
          eeint(ix) = eeint(ix-1) + rmax/(nd-1) * sumx
        enddo
        init=.false.
      endif
      wintsf = polinta(rup,xval,eeint,nd) - polinta(rdn,xval,eeint,nd)
      end

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
c
c when we=w-e ==> 0, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w')
c the integral becomes -v(0)/2
c
c v       = v(iw')
c v0      = v(0)
c x       = s.o.
c wt      = weights for integration
c a       = a constant determining the range of exp(-a^2 w'^2)
c expa(x) = exp(-a^2 w'^2)
c we      = w - e
c nx      = no. x points
      implicit real*8 (a-h,o-z)
      integer:: nx,i,npm
      real(8):: x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx),v0,sum,img=(0d0,1d0),wintz
      real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
      data pi/3.1415926535897932d0/,
cmod 93.06.23
c    d     tol/1.d-6/
     d     tol/1.d-8/
      if (dabs(we) .lt. tol) then
        print *, ' we=',we
        call rx( 'wintz: |w-e| < tol')
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
      integer:: nx,i
      real(8):: x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx),v0,sum
      real(8),parameter :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
      real(8),parameter :: tol=1d-8,pi=4d0*atan(1d0)
      if (dabs(we) .lt. tol) call rx( 'wintz: |w-e| < tol')
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
      if( ratio > 1d0 ) then
        print *, ' ratio =',ratio
        print *, ' gen_ua: ratio >1 strange ratio...'
        call rx( ' gen_ua: ratio >1 strange ratio...')
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
      real(8),intent(in) :: freqx
      real(8),intent(in) :: ratio
      real(8),intent(out) :: ua_
      real(8) :: freqw
      freqw = (1d0 - freqx)/ freqx
      if( ratio > 1d0 ) then
        print *, ' ratio =',ratio
        print *, ' gen_uaa: ratio >1 strange ratio...'
        call rx( ' gen_uaa: ratio >1 strange ratio...')
      else
        ua_ = sqrt(-1d0/freqw*log(ratio))
      endif
      end

!sssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine gen_expa(niw,freqx,ua_,  expa_)
      implicit none
      integer(4), intent(in)  :: niw
      real(8),    intent(in)  :: freqx(niw), ua_
      real(8),    intent(out) :: expa_(niw)
      real(8) :: freqw
      integer(4):: ix
      do ix = 1,niw
        freqw     = (1d0 - freqx(ix))/ freqx(ix)
        expa_(ix) = exp(-(ua_*freqw)**2)
      enddo
      end

!sssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine matzwz2(iSigma_en ,zw,zmel, ntq,nstate,ngb, zwz) !sf 23May02
      use m_mpi, only: mpi__rank
      implicit none
      integer(4), intent(in) :: iSigma_en
      complex(8), intent(in) :: zw(ngb,ngb),zmel(ngb,nstate,ntq)
      integer(4), intent(in) :: nstate,ntq,ngb
      complex(8), intent(out) :: zwz(nstate,ntq,ntq)
      integer(4) :: itp,itpp,it
      complex(8) :: zdotc
      complex(8), allocatable :: CC(:,:,:)
      complex(8),allocatable:: z1r(:,:),z2r(:,:), zwzi(:,:)
      integer :: ivc,verbose
      allocate(CC(ngb,nstate,ntq) )
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntq) !most time consuming part
      if (iSigma_en==1.or.iSigma_en==5) then                                       !sf 23May02
        zwz=0d0
        do itp = 1,ntq
          do  it = 1,nstate
            zwz(it,itp,itp)=
     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itp),1 )
          enddo
        enddo
      elseif (iSigma_en==2 .or. iSigma_en==3 ) then
        if(verbose()>39)write(*,*)'info: USE GEMM FOR SUM (zwz=zmel*cc)'
        allocate(z1r(ntq,ngb),z2r(ngb,ntq),zwzi(ntq,ntq))
        do  it = 1, nstate
          do  itp = 1, ntq
            do  ivc = 1, ngb
              z1r(itp,ivc) = dconjg(zmel(ivc,it,itp))
              z2r(ivc,itp) = CC(ivc,it,itp)
            enddo
          enddo
C         call zprm('z1r',z1r,ntq,ntq,ngb)
C         call zprm('z2r',z2r,ngb,ngb,ntq)
          call zgemm('N','N',ntq,ntq,ngb,(1d0,0d0),z1r,ntq,
     .      z2r,ngb,(0d0,0d0),zwzi,ntq)
          do  itp = 1, ntq
            do itpp = 1, ntq
              zwz(it,itp,itpp) = zwzi(itp,itpp)
            enddo
          enddo
        enddo
        deallocate(z1r,z2r,zwzi)
      else
        call rx( "sxcf, matzwz2: iSigma_en /= 0,1,2,3")
      endif
      deallocate(CC)
      end subroutine matzwz2

      subroutine timeshow(info)
      character*(*) :: info
      write(*,'(a,$)')info
      call cputid(0)
      end
     
!sssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      complex(8) function alagr3zwgt (x,xi,wgt)
      implicit real*8 (a-h,o-z)
c three-point interpolation with arbitrary mesh
c f(x) = [ { (x-x2)(x-x3) } / { (x1-x2)(x1-x3) } ] f1
c      + [ { (x-x1)(x-x3) } / { (x2-x1)(x2-x3) } ] f2
c      + [ { (x-x1)(x-x2) } / { (x3-x1)(x3-x2) } ] f3
c      = wg1(1)* f1+ wgt(2)*f2 +wgt(3)*f3
c x     = the point at which the function is to be interpolated
c xi(3) = points where the function is given
      real(8):: xi(3),wgt(3)
      xx1        = x-xi(1)
      xx2        = x-xi(2)
      xx3        = x-xi(3)
      x12        = xi(1)-xi(2)
      x13        = xi(1)-xi(3)
      x23        = xi(2)-xi(3)
      wgt(1)     =   xx2*xx3/(x12*x13)
      wgt(2)     =   xx1*xx3/(-x12*x23)
      wgt(3)     =   xx1*xx2/(x13*x23)
      end
!> even function version of alagr3z2wgt
      subroutine alagr3z2wgt(x,xi, wgt) 
      implicit none
      intent(in) ::          x,xi
      logical:: ieqj
      real(8) :: amatinv(3,3),amat(3,3),ratio,detxx,wgt(3),x,xi(3)
      if(x<0d0) call rx( ' alagr3z2: x<0d0')
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
      call minv33(amat,amatinv)
      wgt=matmul([1d0,x**2,x**4], amatinv)
      end
      complex(8) function alagr3zz(x,xi,fi) 
! even function version of alagr3z ! return the interpolated value on x for fi(xi).
! Imag part is corrected to be >0
      implicit none
      real(8)::  xi(3), amatinv(3,3),amat(3,3),x,detxx
      complex(8) :: fi(3)
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
      call minv33(amat,amatinv)
      alagr3zz = dcmplx (
     &          sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ),
     &          sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )
      if(dimag(alagr3zz)>0d0) alagr3zz = dcmplx( dreal(alagr3zz),0d0)
      end
