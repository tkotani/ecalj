      subroutine lgen(bas,bmax,nv,nvmax,vecs,work)
c  generates lattice vectors.
      implicit real*8 (a-h,p-z), integer(o)
      dimension bas(3,3),v(3),vecs(3,*),work(*) ! MIZUHO-IR
      call latlim(bas,bmax,imax,jmax,kmax)
      bmax2=bmax*bmax
      nv=0
      do 20 i=-imax,imax
      do 20 j=-jmax,jmax
      do 20 k=-kmax,kmax
        do 21 m=1,3
          v(m)=i*bas(m,1)+j*bas(m,2)+k*bas(m,3)
   21   continue
        v2=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
        if(v2.gt.bmax2) goto 20
        nv=nv+1
        if(nv.gt.nvmax) write(6,633) nvmax,i,imax
        if(nv.gt.nvmax) stop
  633   format(/' --- nv=',i6,'  exceeded,   i=',i3,'  imax=',i3)
        do 22 m=1,3
          vecs(m,nv)=v(m)
   22   continue
        vsm=dabs(v(1))+dabs(v(2))+dabs(v(3))
        work(nv)=v2+vsm/1000.
  20  continue
c --- sort by length -----------
      do 30 iv=1,nv
        ilow=iv
        alow=work(iv)
        do 31 jv=iv,nv
          if(work(jv).lt.alow) then
            alow=work(jv)
            ilow=jv
          endif
  31    continue
        if(ilow.eq.iv) goto 30
        do 32 m=1,3
          xx=vecs(m,iv)
          vecs(m,iv)=vecs(m,ilow)
          vecs(m,ilow)=xx
   32   continue
        work(ilow)=work(iv)
        xx=work(ilow)
c|      write(6,300) iv,(vecs(m,iv),m=1,3),xx
c|300   format(i6,3x,3f9.4,f12.4)
  30  continue
c ---- add neighbor layers if basis vec 3 is not in list ------
      do 41 iv=1,nv
        ddd=(bas(1,3)-vecs(1,iv))**2+(bas(2,3)-vecs(2,iv))**2
     .   +(bas(3,3)-vecs(3,iv))**2
        if(ddd.lt.1.d-8) return
  41  continue
      write(6,650)
  650 format(/' basis vec 3 not in list - include 2 more planes')
      if(3*nv.gt.nvmax) write(6,643) nvmax
      if(3*nv.gt.nvmax) stop
  643 format( '--- lgen needs nvmax at least',i7)
      do 40 iv=1,nv
      do 40 m=1,3
        vecs(m,iv+nv)=vecs(m,iv)+bas(m,3)
        vecs(m,iv+2*nv)=vecs(m,iv)-bas(m,3)
   40 continue
      nv=3*nv
      return
      end
c# latlim fortran *
      subroutine latlim(bas,vmax,i1,i2,i3)
c limits in x y z direction, also initialize matrix
      implicit real*8 (a-h,p-z), integer(o)
      dimension bas(3,3),a(3,3)
      common /clatlm/ a11,a22,a33,a12,a13,a23
      do 6 i=1,3
      do 6 j=i,3
        a(i,j)=bas(1,i)*bas(1,j)+bas(2,i)*bas(2,j)+bas(3,i)*bas(3,j)
    6 continue
      a11=a(1,1)
      a12=a(1,2)
      a13=a(1,3)
      a22=a(2,2)
      a23=a(2,3)
      a33=a(3,3)
      det=a11*a22*a33+a12*a23*a13
     .   +a12*a23*a13-a13*a22*a13
     .   -a23*a23*a11-a12*a12*a33
      i1=vmax*dsqrt((a22*a33-a23**2)/det)
      i2=vmax*dsqrt((a11*a33-a13**2)/det)
      i3=vmax*dsqrt((a11*a22-a12**2)/det)
      return
      end
      subroutine latlm2(vmax,i1,i2,j1,j2,k)
c  input: height in z-direction=k.  output: limits in xy-plane
      implicit real*8 (a-h,p-z), integer(o)
      common /clatlm/ a11,a22,a33,a12,a13,a23
      det2=a11*a22-a12**2
      b11=a22/det2
      b12=-a12/det2
      b22=a11/det2
      b1=a13*k
      b2=a23*k
      xxx=vmax*vmax+b1*b1*b11+b2*b2*b22+2.*b1*b2*b12-a33*k*k
      gx1=-dsqrt(xxx*b11)-(b1*b11+b2*b12)
      gx2= dsqrt(xxx*b11)-(b1*b11+b2*b12)
      gy1=-dsqrt(xxx*b22)-(b1*b12+b2*b22)
      gy2= dsqrt(xxx*b22)-(b1*b12+b2*b22)
      i1=gx1
      if(gx1.gt.0.) i1=i1+1
      i2=gx2
      if(gx2.lt.0.) i2=i2-1
      j1=gy1
      if(gy1.gt.0.) j1=j1+1
      j2=gy2
      if(gy2.lt.0.) j2=j2-1
      return
      end
      subroutine latlm1(vmax,i1,i2,j,k)
c  input: j,k.  output: limits for i.
      implicit real*8 (a-h,p-z), integer(o)
      common /clatlm/ a11,a22,a33,a12,a13,a23
      b=a12*j+a13*k
      c=j*j*a22+k*k*a33+2*j*k*a23-vmax*vmax
      srw=dsqrt(b*b-a11*c)
      g1=(-b-srw)/a11
      g2=(-b+srw)/a11
      i1=g1
      if(g1.gt.0.) i1=i1+1
      i2=g2
      if(g2.lt.0.) i2=i2-1
      return
      end
