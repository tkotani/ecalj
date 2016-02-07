      subroutine lattc(as,tol,alat,alat0,rb0,g1,g2,g3,gt,rb,qb,
     .   lmax,vol,awald,dlat,nkd,rlat,nkr,nkdmx,nkrmx,work)
c  sets up the real and reciprocal space lattice vectors
      implicit real*8 (a-h,p-z), integer(o)
      dimension rlat(3,nkrmx),dlat(3,nkdmx),rb(3,3),qb(3,3),
     .  rb0(3,3),qb0(3,3),work(1) !,ff(0:20),ffb(0:20)
      ipr=iprint()
      tpiba=8.d0*datan(1.d0)/alat
      call cross(rb0(1,2),rb0(1,3),qb0)
      call cross(rb0(1,3),rb0(1,1),qb0(1,2))
      call cross(rb0(1,1),rb0(1,2),qb0(1,3))
      vol0=tripl(rb0,rb0(1,2),rb0(1,3))
      do 34 m=1,3
      do 34 k=1,3
        qb0(m,k)=qb0(m,k)*(1.d0/vol0)
   34 continue
      vol0=dabs(vol0)
      do 11 k=1,3
        call rdist(rb0(1,k),rb(1,k),g1,g2,g3,gt)
        call qdist(qb0(1,k),qb(1,k),g1,g2,g3,gt)
   11 continue
      vol=tripl(rb,rb(1,2),rb(1,3))
      vol=dabs(vol)*(alat**3)

c ... printout
      if (ipr.ge.30) then
        write(6,351)
  351   format(/15x,'real',31x,'recip')
        write(6,350) ((rb0(m,k),m=1,3),(qb0(m,k),m=1,3),k=1,3)
  350   format(3f10.5,5x,3f10.5)
        dstx=dabs(gt-1.d0)
        if(dstx.gt.1d-5) then
          write(6,451) g1,g2,g3,gt
  451     format(/' stretch in direction',3f8.4,'   by',f10.5/)
          write(6,350) ((rb(m,k),m=1,3),(qb(m,k),m=1,3),k=1,3)
          call strain(g1,g2,g3,gt)
        endif
        write(6,998) vol,vol0*(alat**3)
  998   format(' cell volume=',f12.4,6x,'before distortion:',f12.4)
      endif

c ------ set up real and recip vectors ----
      rdist0=vol0**(1.d0/3.d0)
      qdist0=1.d0/rdist0
      radd=.7*rdist0
      qadd=.7*qdist0
      a0=as/rdist0
      awald=a0/alat
      alat1=alat0
      if(alat1.le.0.5d0) alat1=alat
      if(dabs(alat1/alat-1.d0).gt.0.04d0) write(6,560)
  560 format(/' *** warning: alat and alat0 deviate by more than 4 %'/)
      tol1=tol*alat1**(lmax+1)
      call lctoff(a0,vol0,lmax,tol1,r0,q0)
      nkdest=4.18879*(r0+radd)**3/vol0+.5
      nkrest=4.18879*(q0+qadd)**3*vol0+.5
      if (ipr.ge.30) 
     .   write(6,340) as,tol,lmax,awald,vol0,alat1,nkdest,nkrest
  340 format(/' lattc:  as=',f6.3,'   tol=',1p,e8.2,'   lmax=',i1,
     .  '   awald=',0p,f7.4,'   v0=',f10.3/' alat1=',f9.5,
     .  '   estimates:   nkd',i6,'   nkr',i6)
      call lgen(rb0,r0+radd,nkd,nkdmx,dlat,work)
      if (ipr.ge.30) write(6,342) r0,r0*alat,radd,nkd
  342 format('  r0=',f9.4,'   rc=',f9.4,'   radd=',f9.4,'   nkd=', i7)
      call lgen(qb0,q0+qadd,nkr,nkrmx,rlat,work)
      if (ipr.ge.30) write(6,341) q0,q0*tpiba,qadd,nkr
  341 format('  q0=',f9.4,'   qc=',f9.4,'   qadd=',f9.4,'   nkr=', i7)
      do 50 iv=1,nkd
        call rdist(dlat(1,iv),dlat(1,iv),g1,g2,g3,gt)
   50 continue
      do 52 iv=1,nkr
        call qdist(rlat(1,iv),rlat(1,iv),g1,g2,g3,gt)
   52 continue
      return
      end

c --- lctoff
      subroutine lctoff(a0,v0,lmax,tol,r0,q0)
c  makes limits r0,q0 for sums in real and recip space for a lattice
c  with lattice constant 1.
      implicit real*8 (a-h,p-z), integer(o)
      dimension f(0:lmax),g(0:lmax) !kino
      pi=4.d0*datan(1.d0)
      q1=0.001d0
      if(lmax.gt.2) q1=dsqrt(.5d0*(lmax-2))*a0/pi
      gq1=(2d0*pi*q1)**(lmax-2)*dexp(-(pi*q1/a0)**2)*4d0*pi/v0
      if(tol.gt.gq1) write(6,*) '**** lctoff: tol gt gq1'
      q2=50.d0
      q0=5.d0
      do 33 i=1,25
        gq0=(2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
        if(gq0.gt.tol) q1=q0
        if(gq0.lt.tol) q2=q0
        q0=.5d0*(q1+q2)
   33 continue

      r1=0.1d0
      r2=50.d0
      r0=5.d0
      do 15 i=1,25
        call dlmtor(r0,a0,lmax,f,g)
        if(f(lmax).gt.tol) r1=r0
        if(f(lmax).le.tol) r2=r0
        r0=.5d0*(r1+r2)
   15 continue
      try=(2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
c|    write(6,957) q0,try,r0,f(lmax)
c|957 format(' lcut: q0=',f12.6,'   try=',f12.6,'   r0=',f12.6,
c|   .  '   f=',f12.6)
      end

c --- strain
      subroutine strain(gx,gy,gz,gt)
c  printout of strains eps1...eps6 (not engineering strains)
      implicit real*8 (a-h,p-z), integer(o)
      dimension e(3,3),eps(3,3),t(3,3)
      data e/1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/
      do 10 ixyz=1,3
        call rdist(e(1,ixyz),t(1,ixyz),gx,gy,gz,gt)
      do 10 m=1,3
        eps(m,ixyz)=t(m,ixyz)-e(m,ixyz)
   10 continue
      write(6,230) eps(1,1),eps(2,2),eps(3,3),(eps(2,3)+eps(3,2))/2,
     .  (eps(1,3)+eps(3,1))/2,(eps(1,2)+eps(2,1))/2
  230 format(/' tensile strains:     ',3f12.6
     .       /' tensor shear strains:',3f12.6)
      ipr=iprint()
      if(ipr.lt.60) return

c --- print transformation matrix
      det=t(1,1)*t(2,2)*t(3,3)+t(1,2)*t(2,3)*t(3,1)
     .   +t(1,3)*t(2,1)*t(3,2)-t(1,3)*t(2,2)*t(3,1)
     .   -t(1,2)*t(2,1)*t(3,3)-t(1,1)*t(2,3)*t(3,2)
      write(6,441) det,t
  441 format(/' det=',f10.5,'   transf=',3f10.5/(25x,3f10.5))

c --- get expressions for e change to quad order
      alf=gt-1.d0
      if(dabs(alf).lt.1.d-8) return
      f1=eps(1,1)/alf
      f2=eps(2,2)/alf
      f3=eps(3,3)/alf
      f4=0.5d0*(eps(2,3)+eps(3,2))/alf
      f5=0.5d0*(eps(1,3)+eps(3,1))/alf
      f6=0.5d0*(eps(1,2)+eps(2,1))/alf
      write(6,797) f1,f2,f3,f4,f5,f6
  797 format(/' eps/alf= ',6f9.4)

c --- this for cubic crystals
      x11=0.5d0*(f1*f1+f2*f2+f3*f3)
      x12=f1*f2+f1*f3+f2*f3
      x44=0.5d0*(f4*f4+f5*f5+f6*f6)
      write(6,*) 'multipliers for energy change:'
      write(6,944) x11,x12,x44
  944 format(' cubic systems x11,x12,x44',3f9.5)

c --- this for tetragonal systems -------------
      x11=0.5d0*(f1*f1+f2*f2)
      x33=0.5d0*(f3*f3)
      x12=f1*f2
      x13=(f1*f3+f2*f3)
      x44=0.5d0*(f4*f4+f5*f5)
      x66=0.5d0*(f6*f6)
      write(6,945) x11,x33,x12,x13,x44,x66
  945 format(' tetrag: 11,33,12,13,44,66',6f9.5)
      return
      end
