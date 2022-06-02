      subroutine wronkj(e1,e2,r,lmax,fkk,fkj,fjk,fjj)
c  wronskians for hankels and bessels on one sphere.
c  fxy is continuous as e1-->e1.
c  fkk,fjj are symmetric in e1,e2.  fkj(e1,e2)=fjk(e2,e1).
      implicit real*8 (a-h,p-z), integer(o)
      integer:: lmax,lp1,l
      dimension fkk(*),fkj(*),fjk(*),fjj(*),ak1(200),aj1(200), ! MIZUHO-IR
     .   ak2(200),aj2(200),dk2(200),dj2(200),dk1(200),dj1(200)
c ------ first: special case e1=e2=0 -------------
      if(dabs(e1).le.1.d-6.and.dabs(e2).le.1.d-6) then
        r3=r*r*r
        rk=-1.d0
        rj=1.d0/r
        do 20 l=0,lmax
          lp1=l+1
          rk=rk*(2*l-1)/r
          rj=rj*r/(2*l+1)
          fkk(lp1)=rk*rk*r3/(2*l-1)
          fjj(lp1)=-rj*rj*r3/(2*l+3)
          fkj(lp1)=-0.5d0*rj*rk*r3
          fjk(lp1)=fkj(lp1)
   20   continue
        return
      endif
c -------- case e1.ne.e2 --------------------
      if(dabs(e1-e2).gt.1.d-6) then
        efac=1.d0/(e2-e1)
        call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,0)
        call radkj(e2,r,lmax,ak2,aj2,dk2,dj2,0)
        do 10 lp1=1,lmax+1
          fkk(lp1)=efac*r*r*(ak1(lp1)*dk2(lp1)-dk1(lp1)*ak2(lp1))
          fjj(lp1)=efac*r*r*(aj1(lp1)*dj2(lp1)-dj1(lp1)*aj2(lp1))
          fkj(lp1)=efac*r*r*(ak1(lp1)*dj2(lp1)-dk1(lp1)*aj2(lp1))-efac
          fjk(lp1)=efac*r*r*(aj1(lp1)*dk2(lp1)-dj1(lp1)*ak2(lp1))+efac
   10   continue
      else
c --------- case e1.eq.e2 but not zero ---------
        call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,0)
        call radkj(e1,r,lmax,ak2,aj2,dk2,dj2,1)
        do 11 lp1=1,lmax+1
          fkk(lp1)=r*r*(ak1(lp1)*dk2(lp1)-dk1(lp1)*ak2(lp1))
          fjj(lp1)=r*r*(aj1(lp1)*dj2(lp1)-dj1(lp1)*aj2(lp1))
          fkj(lp1)=r*r*(ak1(lp1)*dj2(lp1)-dk1(lp1)*aj2(lp1))
          fjk(lp1)=r*r*(aj1(lp1)*dk2(lp1)-dj1(lp1)*ak2(lp1))
   11   continue
      endif
      return
      end
c --------------------------------------------------------
      subroutine radkj(e,r,lmax,ak,aj,dk,dj,job)
c  radial parts of spherical hankels and bessels.
c  job=0: makes values and slopes. job=1: makes energy derivatives.
      implicit real*8 (a-h,p-z), integer(o)
      integer:: lmax,l,lp1,job
      dimension ak(*),aj(*),dk(*),dj(*),phi(200),psi(200),php(200),psp(200) ! MIZUHO-IR
      er2=e*r*r
      if(job.eq.0) then
        call besslggg(er2,lmax+1,phi,psi)
        rl=1.d0/r
        do 10 l=0,lmax
          lp1=l+1
          rl=rl*r
          ak(lp1)=psi(lp1)/(rl*r)
          aj(lp1)=phi(lp1)*rl
          dk(lp1)=(l*psi(lp1)-psi(l+2))/(rl*r*r)
          dj(lp1)=(l*phi(lp1)-er2*phi(l+2))*rl/r
   10   continue
      else
        call besslggg(er2,lmax+2,phi,psi)
        do 11 lp1=1,lmax+2
          php(lp1)=-0.5d0*r*r*phi(lp1+1)
          psp(lp1)=0.5d0*((2*lp1-1)*psi(lp1)-psi(lp1+1))/e
   11   continue
        rl=1.d0/r
        do 12 l=0,lmax
          lp1=l+1
          rl=rl*r
          ak(lp1)=psp(lp1)/(rl*r)
          aj(lp1)=php(lp1)*rl
          dk(lp1)=(l*psp(lp1)-psp(l+2))/(rl*r*r)
          dj(lp1)=(l*php(lp1)-er2*php(l+2)-r*r*phi(l+2))*rl/r
   12   continue
      endif
      return
      end
