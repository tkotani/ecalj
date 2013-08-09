      subroutine bessl(y,lmax,phi,psi)
c  patched to take care of case lmax=1
c  07.06.91  df floated, to make subroutine work for large lmax's
      implicit real*8 (a-h,p-z), integer(o)
      dimension phi(*),psi(*),dum(420)
      doubleprecision t,dt,t1,tol,my,dt1,dt2
      integer tlp1,tmp1,tmp2
      lmux=max0(lmax,2)
      if(lmux.gt.16) goto 8
      tol=1.d-8
      my=-y
      l=lmux
  1   tlp1=l+l+1
      df=1
      do 2 k=3,tlp1,2
        df=df*k
    2 continue
      t1=1.d0/df
      dt=1.d0
      t=1.d0
      i=0
      do 3 k=1,10000
        i=i+2
        dt1=dble(i)
        dt2=dble(i+tlp1)
        dt=dt*my/(dt1*dt2)
        t=t+dt
        if(dabs(dt).lt.tol) goto 4
  3   continue
      goto 10
  4   if(l.lt.lmux) goto 5
      dum(1)=t1*t
      l=lmux-1
      goto 1
  5   dum(2)=t1*t
      tmp1=lmux+lmux+1
      tmp2=lmux+lmux+2
      nf=tmp1
      do 6 k=3,tmp2
        nf=nf-2
        dum(k)=dble(nf)*dum(k-1)-y*dum(k-2)
    6 continue
      lmuxp1=lmux+1
      lmuxp2=lmux+2
      isn=-1
      do 7 k=1,lmuxp1
        isn=-isn
        j1=lmuxp2-k
        j2=lmuxp1+k
        phi(k)=dum(j1)
        psi(k)=dum(j2)*dble(isn)
    7 continue
      if(y.ge.0.0) return
c ------- negative energy case ----------
      srmy=dsqrt(-y)
      psi(2)=1.d0+srmy
      psi(1)=1.d0
      if(lmux.lt.2) goto 23
      tlp1=1
      do 21 lp1=3,lmuxp1
        tlp1=tlp1+2
        psi(lp1)=dble(tlp1)*psi(lp1-1)-y*psi(lp1-2)
   21 continue
  23  exppr=1.d0/dexp(srmy)
      do 22 lp1=1,lmuxp1
        psi(lp1)=psi(lp1)*exppr
   22 continue
      return
  8   write(6,*) '*** lmux in bessl exceeds max allowed value'
Cstop2rx 2013.08.09 kino      stop
      call rx( '')
  10  write(6,11) y
  11  format(' *** power series in bessl not convergent, y=',e12.4)
Cstop2rx 2013.08.09 kino      stop
      call rx( '')
      end
