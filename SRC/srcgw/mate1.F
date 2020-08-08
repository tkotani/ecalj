c--------------------------------------------------------------------
c mathematics library II
c--------------------------------------------------------------------

      double precision function compvv (v1,v2,n,tol,comment)

c 92.04.18
c compare two vectors v1,v2

      implicit real*8 (a-h,o-z)
      dimension v1(n),v2(n)
      character*72 comment

      do       i = 1,n
        compvv     = dabs(v1(i) - v2(i))
        if (compvv .gt. tol) then
          write (*,6000) comment
Cstop2rx 2013.08.09 kino          stop
          call rx( '')
        endif
      end do

 6000 format (a)
      return
      end
c--------------------------------------------------------------------
      subroutine const (v,c,n)

c 92.03.24
c initialise a vector v of dimension n with a constant c

      implicit real*8 (a-h,o-z)
      dimension v(n)

      do       i = 1,n
        v(i)       = c
      end do

      return
      end
c--------------------------------------------------------------------
      subroutine fbes (z,lmax,
     o                 fs,fc )

c 92.03.06
c generates the coefficients f(n) and f(-n-1) for the
c spherical bessel function for a given argument z
c ref: Abramowitz and Stegun, Handbook of Mathematical Functions
c j(n,z) = f(n,z) sin(z) + (-1)^(n+1) f(-n-1,z) cos(z)
c f(0,z) = 1/z, f(1,z) = 1/z^2

c recurrence relation for f(n,z) is
c f(n-1,z) + f(n+1,z) = f(n,z) (2n+1)/z,  n=0,+/-1,+/-2,....
c f(-n-1,z) + f(-n+1,z) = f(-n,z) (-2n+1)/z

c NOTE:
c problem arises when max(n)/z > 1 since the recurrence starts
c diverging as n ==> very large.
c for small z i.e. large (2n+1)/z, use the small z expansion of j(n,z)
c DO NOT USE THIS ROUTINE

      implicit real*8 (a-h,o-z)
      dimension fs(0:lmax),fc(0:lmax)
      data tol /1.d-75/

Cstop2rx 2013.08.09 kino      if (dabs(z)**(lmax+1) .lt. tol) stop 'fbes: argument z too small'
      if (dabs(z)**(lmax+1) .lt. tol) call rx( 'fbes: argument z too small')
      oz         = 1.d0/z

c sin(z) coefficients, fs(n) = f(n)
      fs(0)      = oz
      fs(1)      = oz*oz
      do       n = 1,lmax-1
        fs(n+1)    = dble(2*n+1) * oz * fs(n) - fs(n-1)
      end do

c cos(z) coefficients, fc(n) = f(-n-1)
c f(-n-1,z) + f(-n+1,z) = f(-n,z) (-2n+1)/z
      fc(0)      = 0.d0
      fc(1)      = -oz
      do       n = 2,lmax
        fc(n)      = dble(1-2*n)*oz*fc(n-1) - fc(n-2)
      end do

c include (-1)^(n+1)
      fsign      = -1.d0
      do       n = 2,lmax
        fc(n)      = fsign*fc(n)
        fsign      = -fsign
      end do

      return
      end
c--------------------------------------------------------------------
      subroutine fbesa (xk,r,
     d                  n,lmax,
     o                  fs,fc )

c 92.03.10
c generate the coefficients of expansions of bessel functions j(l,kr)
c in sin and cos for all r up to lmax, s. fbes

c xk   = k
c r(n) = s.o.
c n    = no. r-points
c lmax = max. l

c fs,fc = expansion coefficents

c NOTE: the coefficients are set to zero if xk*r < tol

      implicit real*8 (a-h,o-z)
      parameter (lmaxim=100)
      dimension r(n)
      dimension fs(n,0:lmax),fc(n,0:lmax)
      dimension wfs(0:lmaxim),wfc(0:lmaxim)
      data tol /1.d-6/

Cstop2rx 2013.08.09 kino      if (lmax .gt. lmaxim) stop 'fbesa: increase lmaxim'
      if (lmax .gt. lmaxim) call rx( 'fbesa: increase lmaxim')

c initialise fs and fc
      do       i = 1,n
        do       l = 0,lmax
          fs(i,l)    = 0.d0
          fc(i,l)    = 0.d0
        end do
      end do

      do 10    i = 1,n
        xkr        = xk*r(i)
        if (dabs(xkr) .lt. tol) goto 10
        call fbes    (xkr,lmax,wfs,wfc)
        do       l = 0,lmax
          fs(i,l)    = wfs(l)
          fc(i,l)    = wfc(l)
        end do
   10 continue

      return
      end
c--------------------------------------------------------------------
      subroutine genbes (lmax,xk,r,
     d                   n,
     o                   bes)

c 92.03.10
c generate spherical bessel functions j(l,kr) up to lmax

c lmax = max. l
c xk   = k, s.o.
c r    = s.o.
c n    = no. r-points

c bes(r,l) = spherical bessel functions

      implicit real*8 (a-h,o-z)
      parameter (lmaxim=100)
      dimension r(n)
      dimension bes(n,0:lmax)
      dimension spb(0:lmaxim)

Cstop2rx 2013.08.09 kino      if (lmax .gt. lmaxim) stop 'genbes: increase lmaxim'
      if (lmax .gt. lmaxim) call rx( 'genbes: increase lmaxim')

      xk         = dabs(xk)
      do       i = 1,n
        rk         = xk*r(i)
        call sphbes (rk,lmax,spb)
        do       l = 0,lmax
          bes(i,l)   = spb(l)
        end do
      end do

      return
      end
c--------------------------------------------------------------------
      integer function icompvv (iv1,iv2,n,comment)

c 92.04.18
c compare two integer vectors iv1,iv2

      implicit real*8 (a-h,o-z)
      dimension iv1(n),iv2(n)
      character*72 comment

      do       i = 1,n
        icompvv    = iv1(i) - iv2(i)
        if (icompvv .ne. 0) then
          write (*,6000) comment
Cstop2rx 2013.08.09 kino          stop
          call rx( '')
        endif
      end do

 6000 format (a)
      return
      end
c--------------------------------------------------------------------
      subroutine mtrans (a,ldim1,ldim2,at)

c 93.04.12 transpose a matrix at(i,j) = a(j,i)

      implicit real*8 (a-h,o-z)
      dimension a(ldim1,ldim2),
     o          at(ldim2,ldim1)

      do       j = 1,ldim1
        do       i = 1,ldim2
          at(i,j)    = a(j,i)
        enddo
      enddo
      return
      end
c--------------------------------------------------------------------
      subroutine mvc (rm,cm,rv,cv,code,
     d                ldm,lrow,nsum,
     o                rmv,cmv)

c 92.08.20
c multiply a complex matrix m and a complex vector v
c if code = 'N' or other than 'C':
c mv(i) = S[k=1,nsum] m(i,k) v(k),  i=1,lrow
c if code = 'C':
c mv(i) = S[k=1,nsum] m(i,k) cc[v(k)],  i=1,lrow

c rm,cm = Re and Im matrix m
c rv,cv = Re and Im vector v
c ldm   = leading dimension of m
c lrow  = s. above
c nsum  = s. above

c rmv,cmv = Re and Im mv

      implicit real*8 (a-h,o-z)
      character*1 code
      dimension rm(ldm,nsum),cm(ldm,nsum),
     i          rv(nsum),cv(nsum),
     o          rmv(lrow),cmv(lrow)

      if (code .eq. 'C') goto 1111
      do       i = 1,lrow
        rsum       = 0.d0
        csum       = 0.d0
        do       k = 1,nsum
          rsum       = rsum + rm(i,k)*rv(k) - cm(i,k)*cv(k)
          csum       = csum + rm(i,k)*cv(k) + cm(i,k)*rv(k)
        enddo
        rmv(i)     = rsum
        cmv(i)     = csum
      enddo
      return

 1111 do       i = 1,lrow
        rsum       = 0.d0
        csum       = 0.d0
        do       k = 1,nsum
          rsum       = rsum + rm(i,k)*rv(k) + cm(i,k)*cv(k)
          csum       = csum - rm(i,k)*cv(k) + cm(i,k)*rv(k)
        enddo
        rmv(i)     = rsum
        cmv(i)     = csum
      enddo

      return
      end
c--------------------------------------------------------------------
      double precision function simpsx (f,a,nr)
      implicit real*8(a-h,o-z)

c integrate dr f(r) with equal mesh a, nr points including end points,
c with Simpson's rule.

      dimension f(nr)
      data third/0.3333333333333333d0/

cmod 1998.09.22
      simpsx     = 0.d0
      if (nr .le. 1) return
      if (nr .eq. 2) then
        simpsx     = 0.5d0*a*(f(1)+f(2))
        return
      endif

      nr2        = nr/2
      nrt        = nr
      if(nr2*2 .eq. nr)nrt = nr - 1
      nr2        = nrt - 1

      even       = 0.d0
      odd        = 0.d0
      do 10    i = 2,nr2,2
        even       = even + f(i)
        odd        = odd  + f(i+1)
   10 continue
      result     = f(1) - f(nrt)
     +           + 4.d0*even + 2.d0*odd

      if(nrt .lt. nr)
     +result     = result
     +    + 0.25d0*(      -f(nr-2)
     +               +8.d0*f(nr-1)
     +               +5.d0*f(nr)  )

      simpsx     = result*a*third

      return
      end
c---------------------------------------------------------------------
      double precision function simpsr (f,a,br,nr)

c 92.03.11
c integrate dr f(r) with r(i) = b*( exp((i-1)*a) - 1 )
c br = b + r
c with Simpson's rule.

      implicit real*8(a-h,o-z)
      dimension f(nr),br(nr)
      data third/0.3333333333333333d0/

      nr2        = nr/2
      nrt        = nr
      if(nr2*2 .eq. nr)nrt = nr - 1
      nr2        = nrt - 1

      even       = 0.d0
      odd        = 0.d0
      do 10    i = 2,nr2,2
        even       = even + br(i) * f(i)
        odd        = odd  + br(i+1) * f(i+1)
   10 continue
      result     = br(1)*f(1) - br(nrt)*f(nrt)
     +           + 4.d0*even + 2.d0*odd

      if(nrt .lt. nr)
     +result     = result
     +    + 0.25d0*(      -br(nr-2)*f(nr-2)
     +               +8.d0*br(nr-1)*f(nr-1)
     +               +5.d0*br(nr)  *f(nr)  )

      simpsr     = result*a*third
      return
      end
c---------------------------------------------------------------------
      double precision function simpsr2(f,g,a,br,nr)

c 92.08.24 fr sub.simpson2
c integrate dr f(r)*g(r) with r(i) = b*( exp((i-1)*a) - 1 )
c and br = b+r
c with Simpson's rule.

      implicit real*8(a-h,o-z)
      dimension f(nr),g(nr),br(nr)
      data third/0.3333333333333333d0/

      nr2        = nr/2
      nrt        = nr
      if(nr2*2 .eq. nr)nrt = nr - 1
      nr2        = nrt - 1

      even       = 0.d0
      odd        = 0.d0
      do 10    i = 2,nr2,2
        even       = even + br(i) * f(i) * g(i)
        odd        = odd  + br(i+1) * f(i+1) * g(i+1)
   10 continue
      result     = br(1)*f(1)*g(1) 
     +           - br(nrt)*f(nrt)*g(nrt)
     +           + 4.d0*even + 2.d0*odd

      if(nrt .lt. nr)
     +result     = result
     +    + 0.25d0*(      -br(nr-2)*f(nr-2)*g(nr-2)
     +               +8.d0*br(nr-1)*f(nr-1)*g(nr-1)
     +               +5.d0*br(nr)  *f(nr)  *g(nr)  )

      simpsr2    = result*a*third
      return
      end
c---------------------------------------------------------------------
      double precision function vdot4 (f1,f2,f3,f4,n)

c calculates vdot4 = S[i=1,n] f1(i)*f2(i)*f3(i)*f4(i)

      implicit real*8 (a-h,o-z)
      dimension f1(n),f2(n),f3(n),f4(n)

      vdot4      = 0.d0
      do       i = 1,n
        vdot4      = vdot4 + f1(i)*f2(i)*f3(i)*f4(i)
      end do

      return
      end
c--------------------------------------------------------------------
      double precision function vmax (v,n)

c 92.05.18
c finds the maximum element in a vector v

      implicit real*8 (a-h,o-z)
      dimension v(n)

      vmax       = v(1)
      do       i = 2,n
        if (v(i) .gt. vmax) vmax = v(i)
      enddo

      return
      end
c--------------------------------------------------------------------
      double precision function xlorenz (w,delta)

c 93.05.19 lorentzian = (1/pi) delta/[ w^2 + delta^2 ]
c normalized to one

      implicit real*8 (a-h,o-z)
      data pi /3.1415926535897932/
      xlorenz    = delta/(pi*(w*w+delta*delta))
      return
      end
c--------------------------------------------------------------------
