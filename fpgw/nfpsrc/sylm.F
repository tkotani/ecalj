      subroutine sylm(r,yl,lmx,r2s)
c  unnormalized spherical harmonic polynomials (from asw progs)
c  polar axis along 001 axis. all in double precision
      implicit real*8 (a-h,p-z), integer(o)
      dimension r(3),yl(*) !,c(15),s(15),p(15,15) !kino yl(4) --->yl(*)
      integer(4),parameter:: lmxx=50
      double precision c(lmxx+1),s(lmxx+1),p(lmxx+1,lmxx+1)
      equivalence (x,c(2)),(y,s(2)),(z,p(2,1))
      data c(1),s(1),p(1,1),p(2,2) /1.d0,0.d0,1.d0,1.d0/
      n=(lmx+1)**2
      yl(1)=1.d0
      x=r(1)
      y=r(2)
      z=r(3)
  500 st=x*x+y*y
      z2=z*z
      r2=st+z2
      r2s=r2
      if(n.lt.2) return
      if(r2.gt.1.d-10) goto 1
      do 6 i=2,n
        yl(i)=0.d0
    6 continue
      return
  1   yl(2)=y
      yl(3)=z
      yl(4)=x
      nt=1
      do 5 l=2,lmx
        lp1=l+1
        lm1=l-1
        lav=l*lp1+1
        p(lp1,1)=(dfloat(l+lm1)*z*p(l,1)-dfloat(lm1)*r2*p(lm1,1))
     .   /dfloat(l)
        yl(lav)=p(lp1,1)
        nt=nt+2
        p(lp1,lp1)=p(l,l)*dfloat(nt)
        c(lp1)=x*c(l)-y*s(l)
        s(lp1)=x*s(l)+y*c(l)
        lavpl=lav+l
        yl(lavpl)=p(lp1,lp1)*c(lp1)
        lavml=lav-l
        yl(lavml)=p(lp1,lp1)*s(lp1)
        if(st.gt.z2) goto 3
        do 2 lmm=1,lm1
          m=l-lmm
          lavpm=lav+m
          lavmm=lav-m
          mp1=m+1
          p(lp1,mp1)=(r2*dfloat(l+m)*p(l,mp1)-st*p(lp1,mp1+1))
     .  /(z*dfloat(l-m))
          yl(lavpm)=p(lp1,mp1)*c(mp1)
          yl(lavmm)=p(lp1,mp1)*s(mp1)
    2   continue
        goto 5
  3     do 4 m=1,lm1
          mp1=m+1
          lavpm=lav+m
          lavmm=lav-m
          p(lp1,mp1)=(dfloat(lm1+m)*r2*p(l,m)-dfloat(lp1-m)*z*p(lp1,m))/st
          yl(lavpm)=p(lp1,mp1)*c(mp1)
  4       yl(lavmm)=p(lp1,mp1)*s(mp1)
  5     continue
        return
        end
