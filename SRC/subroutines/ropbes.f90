module m_ropbes !  Spherical Bessel function
  public ropbes
  private
  contains
! ----------------------------------------------------------------
!   spherical Bessel function at x=r(i)*sqrt(e) divided by x**l
!i Inputs
!i   r    list of points
!i   e    energy
!i   y,h  work vectors of length n each (not used)
!o Outputs
!o   xi   J(r,l)/r**l, according to standard definition
! origianl: Feb. 15, 2010, Hiori Kino
subroutine ropbes(r,e,lmax,y,h,xi,n)
  implicit none
  integer :: lmax,n
  double precision :: e,r(n),xi(n,0:lmax),h(n),y(n)
  double precision :: e2,f(0:lmax),eps=1.0d-10,x
  double precision :: sj,sy,sjp,syp
  integer:: i,l
  e2=sqrt(e)
  do i=1,n
     x=r(i)*e2
     if (x< eps)  x=eps
     do l=0,lmax
        call sphbes(l,x,sj,sy,sjp,syp)
        f(l)=sj/(x**l)
     enddo
     xi(i,0:lmax)= f(0:lmax)
  enddo
end subroutine ropbes
!--------------------------------
SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
  implicit none
  INTEGER :: MAXIT
  REAL(8) rj,rjp,ry,ryp,x,xnu,XMIN
  DOUBLE PRECISION :: EPS,FPMIN,PI
  PARAMETER (EPS=1.e-10,FPMIN=1.e-30,MAXIT=10000,XMIN=2., &
       PI=3.141592653589793d0)
  !U    USES beschb
  INTEGER :: i,isign,l,nl
  DOUBLE PRECISION :: a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e, &
       f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q, &
       r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1, &
       temp,w,x2,xi,xi2,xmu,xmu2
  if(x <= 0. .OR. xnu < 0.) then
     write(*,*) 'bad arguments in bessjy',x,xnu
     stop
  endif
  if(x < XMIN)then
     nl=int(xnu+.5d0)
  else
     nl=max(0,int(xnu-x+1.5d0))
  endif
  xmu=xnu-nl
  xmu2=xmu*xmu
  xi=1.d0/x
  xi2=2.d0*xi
  w=xi2/PI
  isign=1
  h=xnu*xi
  if(h < FPMIN)h=FPMIN
  b=xi2*xnu
  d=0.d0
  c=h
  do 11 i=1,MAXIT
     b=b+xi2
     d=b-d
     if(abs(d) < FPMIN)d=FPMIN
     c=b-1.d0/c
     if(abs(c) < FPMIN)c=FPMIN
     d=1.d0/d
     del=c*d
     h=del*h
     if(d < 0.d0)isign=-isign
     if(abs(del-1.d0) < EPS)goto 1
11 enddo
  write(*,*) 'x too large in bessjy; try asymptotic expansion'
  stop
1 continue
  rjl=isign*FPMIN
  rjpl=h*rjl
  rjl1=rjl
  rjp1=rjpl
  fact=xnu*xi
  do 12 l=nl,1,-1
     rjtemp=fact*rjl+rjpl
     fact=fact-xi
     rjpl=fact*rjtemp-rjl
     rjl=rjtemp
12 enddo
  if(rjl == 0.d0)rjl=EPS
  f=rjpl/rjl
  if(x < XMIN) then
     x2=.5d0*x
     pimu=PI*xmu
     if(abs(pimu) < EPS)then
        fact=1.d0
     else
        fact=pimu/sin(pimu)
     endif
     d=-log(x2)
     e=xmu*d
     if(abs(e) < EPS)then
        fact2=1.d0
     else
        fact2=sinh(e)/e
     endif
     call beschb(xmu,gam1,gam2,gampl,gammi)
     ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
     e=exp(e)
     p=e/(gampl*PI)
     q=1.d0/(e*PI*gammi)
     pimu2=0.5d0*pimu
     if(abs(pimu2) < EPS)then
        fact3=1.d0
     else
        fact3=sin(pimu2)/pimu2
     endif
     r=PI*pimu2*fact3*fact3
     c=1.d0
     d=-x2*x2
     sum=ff+r*q
     sum1=p
     do 13 i=1,MAXIT
        ff=(i*ff+p+q)/(i*i-xmu2)
        c=c*d/i
        p=p/(i-xmu)
        q=q/(i+xmu)
        del=c*(ff+r*q)
        sum=sum+del
        del1=c*p-i*del
        sum1=sum1+del1
        if(abs(del) < (1.d0+abs(sum))*EPS)goto 2
13   enddo
     write(*,*) 'bessy series failed to converge'
     stop
2    continue
     rymu=-sum
     ry1=-sum1*xi2
     rymup=xmu*xi*rymu-ry1
     rjmu=w/(rymup-f*rymu)
  else
     a=.25d0-xmu2
     p=-.5d0*xi
     q=1.d0
     br=2.d0*x
     bi=2.d0
     fact=a*xi/(p*p+q*q)
     cr=br+q*fact
     ci=bi+p*fact
     den=br*br+bi*bi
     dr=br/den
     di=-bi/den
     dlr=cr*dr-ci*di
     dli=cr*di+ci*dr
     temp=p*dlr-q*dli
     q=p*dli+q*dlr
     p=temp
     do 14 i=2,MAXIT
        a=a+2*(i-1)
        bi=bi+2.d0
        dr=a*dr+br
        di=a*di+bi
        if(abs(dr)+abs(di) < FPMIN)dr=FPMIN
        fact=a/(cr*cr+ci*ci)
        cr=br+cr*fact
        ci=bi-ci*fact
        if(abs(cr)+abs(ci) < FPMIN)cr=FPMIN
        den=dr*dr+di*di
        dr=dr/den
        di=-di/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        if(abs(dlr-1.d0)+abs(dli) < EPS)goto 3
14   enddo
     write(*,*) 'cf2 failed in bessjy'
     stop
3    continue
     gam=(p-f)/q
     rjmu=sqrt(w/((p-f)*gam+q))
     rjmu=sign(rjmu,rjl)
     rymu=rjmu*gam
     rymup=rymu*(p+q/gam)
     ry1=xmu*xi*rymu-rymup
  endif
  fact=rjmu/rjl
  rj=rjl1*fact
  rjp=rjp1*fact
  do 15 i=1,nl
     rytemp=(xmu+i)*xi2*ry1-rymu
     rymu=ry1
     ry1=rytemp
15 enddo
  ry=rymu
  ryp=xnu*xi*rymu-ry1
  return
end SUBROUTINE bessjy
!  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
  implicit none
  INTEGER :: n
  REAL(8) sj,sjp,sy,syp,x
  !U    USES bessjy
  REAL(8) factor,order,rj,rjp,ry,ryp,RTPIO2
  PARAMETER (RTPIO2=1.2533141)
  if(n < 0 .OR. x <= 0.)then
     write(*,*)  'bad arguments in sphbes'
     stop
  endif
  order=n+0.5
  call bessjy(x,order,rj,ry,rjp,ryp)
  factor=RTPIO2/sqrt(x)
  sj=factor*rj
  sy=factor*ry
  sjp=factor*rjp-sj/(2.*x)
  syp=factor*ryp-sy/(2.*x)
  return
end SUBROUTINE sphbes
SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
  implicit none
  INTEGER :: NUSE1,NUSE2
  DOUBLE PRECISION :: gam1,gam2,gammi,gampl,x
  PARAMETER (NUSE1=5,NUSE2=5)
  REAL(8) :: xx,c1(7),c2(8)
  SAVE c1,c2
  DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4, &
       -3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
  DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3, &
       -4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
  xx=8.d0*x*x-1.d0
  gam1=chebev(-1d0,1d0,c1,NUSE1,xx)
  gam2=chebev(-1d0,1d0,c2,NUSE2,xx)
  gampl=gam2-x*gam1
  gammi=gam2+x*gam1
  return
end SUBROUTINE beschb
!  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
real(8) FUNCTION chebev(a,b,c,m,x)
  INTEGER :: m
  REAL(8) :: a,b,x,c(m)
  INTEGER :: j
  REAL(8) :: d,dd,sv,y,y2 !real(8) 2022-6-6
  if ((x-a)*(x-b) > 0d0)then
     write(*,*) 'x not in range in chebev'
     stop
  endif
  d=0d0
  dd=0d0
  y=(2d0*x-a-b)/(b-a)
  y2=2d0*y
  do 11 j=m,2,-1
     sv=d
     d=y2*d-dd+c(j)
     dd=sv
11 enddo
  chebev=y*d-dd+.5d0*c(1)
  return
end FUNCTION chebev

end module m_ropbes
