      subroutine rdeq(enum1,enum2,ksop,z,v,rofi,nr,nsp,a,b,l,lmx,
     .                imu,scalede,gmt,fmt,gmtde,fmtde,gsmt,pprel)

C-Solves radial dirac equations, in the presence of an effective magnetic field
C ----------------------------------------------------------------------
Ci Inputs
Ci   enum1 :?
Ci   enum2 :?
Ci   ksop  :?
Ci   z     :nuclear charge
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   l     :l quantum number
Ci   lmx   :dimensions ksop and pprel
Ci   imu   :
Ci   scalede:
Ci   gmt   :
Ci   fmt   :
Ci   gsmt  :
Ci   pprel :
Co Outputs
Co   gmtde :(2x2 matrix) r*g (g=normalized large component) at MT sphere
Co   fmtde :(2x2 matrix) r*f (f=normalized small component) at MT sphere
Cl Local variables
Cr Remarks
Cr   In this version enu depends only on kappa
Cu Updates
Cu   18 Jun 04 (A Chantis) working version
Cu   24 Jun 03 (A Chantis) adopted from V. Antropov
C ----------------------------------------------------------------------
C     implicit none
C     Passed parameters
      integer  nr, nsp, imu, l, lmx
      double precision v(nr,nsp), enu, rofi(nr), z,
     .  a, b, enum1, enum2, scalede
C     Local variables
      double precision xx1(4), xx2(4), dx, hoc, cc, xk, xk1, enu1, enu2,
     .                 nug1(nr), nug2(nr), nuf1(nr), nuf2(nr),
     .                 dg11, df11, dg21, df21, dg12, df12, dg22,
     .                 df22, dg13, df13, dg23, df23, dg14, df14,
     .                 dg24, df24, u, up, upp, J, q, r, bm(nr),
     .                 vm(nr), e1,e2,vmint, bmint, g1c, f1c, g2c,
     .                 f2c, nug1int, nuf1int, nug2int,
     .                 nuf2int, de, gamma, csq, g1sn, f1sn,
     .                 g1snm1, g1snm2, g1snm3, g1snm4, f1snm1,
     .                 f1snm2, f1snm3, f1snm4, g2sn, f2sn, g2snm1,
     .                 g2snm2, g2snm3, g2snm4, f2snm1, f2snm2,
     .                 f2snm3, f2snm4, g1sng1, f1snf1, g2sng2, f2snf2
      double precision 
     .                 g1p, f1p, g2p, f2p, norm, gmt(2,2),
     .                 fmt(2,2), gmt0(2,2,5), fmt0(2,2,5), gmtde(2,2),
     .                 fmtde(2,2),smalpa(2,2),gn1(2,5,nr),gn2(2,5,nr),
     .                 gn1dot(2,nr), gn2dot(2,nr),fn1(2,5,nr),
     .                 fn2(2,5,nr), fn1dot(2,nr), fn2dot(2,nr),
     .                 factor(nr), factor1(nr),  factor2(nr),
     .                 factor3(nr), product, product1,
     .                 gsmt(2,2), psi(4,nr), ksi(4,nr), mu,
     .                 ksop(0:lmx,nsp,nsp,6), smalpa1(2,2)
C     double precision factor4(nr),factor5(nr)
      integer n, ie, alpha, alp1, alp2, nrk, ir,
     .        nit, nrmx, i, k
C     Output
      double precision g1(nr), g2(nr), f1(nr), f2(nr),
     .                 g1s(nr), g2s(nr), f1s(nr), f2s(nr)
      double precision pprel(4,0:lmx,2*(lmx+1),2,2)

C     open(78,file='test.dat',status='unknown')

cc      print*, 'THE B and NR RDEQ gets', b, nr, mu
c      if (l .eq. 2) then
      enu1 = -0.5d0*dble((l+1))*ksop(l,1,1,1) + 0.5d0*(enum1 + enum2)
      enu2 = 0.5d0*dble(l)*ksop(l,1,1,1) + 0.5d0*(enum1 + enum2)
c      else
c      enu1 = 0.5d0*(enum1 + enum2)
c      enu2 = 0.5d0*(enum1 + enum2)
c      endif
c      print*,'====================='
c      print*, l,ksop(l,1,1,1),enum1,enum2
c      print*, l,enu1,enu2

      enu = 0.d0

      mu = dble(imu-l) - 1.5d0

cc---Test: Eliminate Spin-Orbit for non-mag. materials
cc     mu = dble(l)+0.5d0

cc	enu = enum1
cc      write(78,*) l,enum1
cc      goto 13

      if (z .eq. 0) return

      cc = 274.072d0
      csq = cc*cc
      de = 0.03d0*scalede
cc      de = 5.d-10
      nrk = 6
      nrmx = nr
      hoc = (2.d0*z/cc)
      dx = 1.0d0

cc      do IAT= 1, NTYPE

cc      do l= 0, lmax

cc      mumax= 2*(l+1)
cc      print *, mumax

cc      do imu= 1, mumax
cc        imu = 1
cc         mu= dble(imu - l) - 1.5d0

      if (dabs(mu) .eq. (l+0.5d0)) then
        xk1= dble(-l-1)
        xk = xk1
      else
        xk1=dble(-l-1)
        xk =dble(l)
      endif

C----Initial value: derived from initial conditions--
C     r-->0, V = V0 - 2Z/r, B = const
C     The wave function is a polymonial
C     g0i = r^(gi-1)*Sum_over_v(g1v*r^v)
C----------------------------------------------------
      gamma = dsqrt(xk*xk - hoc*hoc)
      xx1(1) = 1.d0*(rofi(2)**(gamma-1.d0))
      xx1(2) = (xk + gamma)/hoc*(rofi(2)**(gamma-1.d0))
      xx1(3) = 0.d0
      xx1(4) = 0.d0

      gamma = dsqrt(xk1*xk1 - hoc*hoc)
      xx2(1) = 0.d0
      xx2(2) = 0.d0
      xx2(3) = 1.d0*(rofi(2)**(gamma-1.d0))
      xx2(4) = (xk1 + gamma)/hoc*(rofi(2)**(gamma-1.d0))

      do ir = 2, nr
        vm(ir) = 0.5d0*(v(ir,1) + v(ir,2))
        bm(ir) = 0.5d0*(v(ir,2) - v(ir,1))
cc        vm(ir) = v(ir,1)
cc        bm(ir) = 0.0d0
      enddo

      up = mu/(l-0.5d0)
      u = mu/(l+0.5d0)
      upp = mu/(l+1.5d0)


C-----Added to take care of the abs(mu)=l+0.5 case------
      if(dabs(mu) .eq. (l+0.5d0)) then
        up = upp
      endif
C-------------------------------------------------------

      g1(1) = 0.0d0
      f1(1) = 0.0d0
      g2(1) = 0.0d0
      f2(1) = 0.0d0

      do alpha = 1, 2
        if (alpha .eq. 2) go to 100
        g1(2) = xx1(1)
        f1(2) = xx1(2)
        g2(2) = xx1(3)
        f2(2) = xx1(4)
        go to 102
 100    g1(2) = xx2(1)
        f1(2) = xx2(2)
        g2(2) = xx2(3)
        f2(2) = xx2(4)

 102    do ie= 1, 5
          e2= enu2 + (ie - 3)*de
          e1= enu1 + (ie - 3)*de
cc 102     e = enu
cc         print *, 'e-v=', e-vm(2)

          do ir = 2, nr
            nug1(ir) = 1.d0 + (1.d0/csq*(e1 + 2.d0*z/rofi(ir)
     .            - vm(ir) + up*bm(ir)))

            nuf1(ir) = - (e1 + 2.d0*z/rofi(ir) - vm(ir) + u*bm(ir))

C---Added to take care of the abs(mu)=l+0.5 case-----------
            if(dabs(mu) .eq. (l+0.5d0)) then
              nug1(ir) = 1.d0 + (1.d0/csq*(e2 + 2.d0*z/rofi(ir)
     .            - vm(ir) - upp*bm(ir)))

              nuf1(ir) = - (e2 + 2.d0*z/rofi(ir) - vm(ir) - u*bm(ir))
            endif
C----------------------------------------------------------

            nug2(ir) = 1.d0 + (1.d0/csq*(e2 + 2.d0*z/rofi(ir)
     .            - vm(ir) - upp*bm(ir)))

            nuf2(ir) = - (e2 + 2.d0*z/rofi(ir) - vm(ir) - u*bm(ir))
          enddo

          J = a*(rofi(2) + b)
          q = dexp(a/2.0d0)
          r = rofi(2)

          g1s(2) = J*(cc*nug1(2)*f1(2) - xk/r*g1(2))
          f1s(2) = J*((nuf1(2)*g1(2) - dsqrt(1.d0 - u**2)*bm(2)*g2(2))/cc
     .  + xk/r*f1(2))
          g2s(2) = J*(cc*nug2(2)*f2(2) - xk1/r*g2(2))
          f2s(2) = J*((nuf2(2)*g2(2) - dsqrt(1.d0 - u**2)*bm(2)*g1(2))/cc
     .  + xk1/r*f2(2))

C-----------Runge-Kutta for the first NKR points ------------------------------


C-----First point--------------------------------------------------------------

          n = 2

 2        g1c = g1(n)
          f1c = f1(n)
          g2c = g2(n)
          f2c = f2(n)

          dg11 = dx*J*(cc*nug1(n)*f1c - xk/r*g1c)
          df11 = dx*J*((nuf1(n)*g1c - dsqrt(1.d0 - u**2)*bm(n)*g2c)/cc
     .  + xk/r*f1c)
          dg21 = dx*J*(cc*nug2(n)*f2c - xk1/r*g2c)
          df21 = dx*J*((nuf2(n)*g2c - dsqrt(1.d0 - u**2)*bm(n)*g1c)/cc
     .  + xk1/r*f2c)

          g1c = g1(n) + 0.5d0*dg11
          f1c = f1(n) + 0.5d0*df11
          g2c = g2(n) + 0.5d0*dg21
          f2c = f2(n) + 0.5d0*df21

          J = J*q
          r = J/a - b

          vmint = (5.d0*vm(n) + 2.d0*vm(n+1) + vm(n+2))*.125d0
          bmint = (5.d0*bm(n) + 2.d0*bm(n+1) + bm(n+2))*.125d0

cc      vmint = (6.d0*vm(n) + 5.d0*vm(n+1) - 3.d0*vm(n-1))*.125d0
cc      bmint = (6.d0*bm(n) + 5.d0*bm(n+1) - 3.d0*vm(n-1))*.125d0


          nug1int = 1.d0 + (1.d0/csq*(e1  + 2.d0*z/r
     .          - vmint + up*bmint))

          nuf1int = - (e1 + 2.d0*z/r - vmint + u*bmint)

C---Added to take care of the abs(mu)=l+0.5 case-----------
          if(dabs(mu) .eq. (l+0.5d0)) then
            nug1int = 1.d0 + (1.d0/csq*(e2 + 2.d0*z/r
     .          - vmint - upp*bmint))
            nuf1int = - (e2 + 2.d0*z/r - vmint - u*bmint)
          endif
C---------------------------------------------------------

          nug2int = 1.d0 + (1.d0/csq*(e2 + 2.d0*z/r
     .          - vmint - upp*bmint))
          nuf2int = - (e2 + 2.d0*z/r - vmint - u*bmint)


          dg12 = dx*J*(cc*nug1int*f1c - xk/r*g1c)
          df12 = dx*J*((nuf1int*g1c - dsqrt(1.d0 - u**2)*bmint*g2c)/cc
     .  + xk/r*f1c)
          dg22 = dx*J*(cc*nug2int*f2c - xk1/r*g2c)
          df22 = dx*J*((nuf2int*g2c - dsqrt(1.d0 - u**2)*bmint*g1c)/cc
     .  + xk1/r*f2c)

          g1c = g1(n) + 0.5d0*dg12
          f1c = f1(n) + 0.5d0*df12
          g2c = g2(n) + 0.5d0*dg22
          f2c = f2(n) + 0.5d0*df22

          dg13 = dx*J*(cc*nug1int*f1c - xk/r*g1c)
          df13 = dx*J*((nuf1int*g1c - dsqrt(1.d0 - u**2)*bmint*g2c)/cc
     .  + xk/r*f1c)
          dg23 = dx*J*(cc*nug2int*f2c - xk1/r*g2c)
          df23 = dx*J*((nuf2int*g2c - dsqrt(1.d0 - u**2)*bmint*g1c)/cc
     .  + xk1/r*f2c)

          g1c = g1(n) + dg13
          f1c = f1(n) + df13
          g2c = g2(n) + dg23
          f2c = f2(n) + df23

          n = n + 1

          J = J*q
          r = J/a - b

          dg14 = dx*J*(cc*nug1(n)*f1c - xk/r*g1c)
          df14 = dx*J*((nuf1(n)*g1c - dsqrt(1.d0 - u**2)*bm(n)*g2c)/cc
     .  + xk/r*f1c)
          dg24 = dx*J*(cc*nug2(n)*f2c - xk1/r*g2c)
          df24 = dx*J*((nuf2(n)*g2c - dsqrt(1.d0 - u**2)*bm(n)*g1c)/cc
     .  + xk1/r*f2c)

          g1(n) = g1(n-1) + (dg11 + 2.d0*(dg12+dg13) + dg14)/6.d0
          f1(n) = f1(n-1) + (df11 + 2.d0*(df12+df13) + df14)/6.d0
          g2(n) = g2(n-1) + (dg21 + 2.d0*(dg22+dg23) + dg24)/6.d0
          f2(n) = f2(n-1) + (df21 + 2.d0*(df22+df23) + df24)/6.d0

C----Determine the derivatives dP/dx and dQ/dx----------------------------

          g1s(n) = J*(cc*nug1(n)*f1(n) - xk/r*g1(n))
          f1s(n) = J*((nuf1(n)*g1(n) - dsqrt(1.d0 - u**2)*bm(n)*g2(n))/cc
     .  + xk/r*f1(n))
          g2s(n) = J*(cc*nug2(n)*f2(n) - xk1/r*g2(n))
          f2s(n) = J*((nuf2(n)*g2(n) - dsqrt(1.d0 - u**2)*bm(n)*g1(n))/cc
     .  + xk1/r*f2(n))

C---------------------------------------------------------------------
C     The rest of the integration is performed using Milne's method.
C---------------------------------------------------------------------

          if(n-nrk) 2,3,3

 3        if(n .eq. nrmx) goto 12

          g1sn = g1s(n)
          f1sn = f1s(n)
          g1snm1 = g1s(n-1)
          f1snm1 = f1s(n-1)
          g1snm2 = g1s(n-2)
          f1snm2 = f1s(n-2)
          g1snm3 = g1s(n-3)
          f1snm3 = f1s(n-3)
          g1snm4 = g1s(n-4)
          f1snm4 = f1s(n-4)

          g2sn = g2s(n)
          f2sn = f2s(n)
          g2snm1 = g2s(n-1)
          f2snm1 = f2s(n-1)
          g2snm2 = g2s(n-2)
          f2snm2 = f2s(n-2)
          g2snm3 = g2s(n-3)
          f2snm3 = f2s(n-3)
          g2snm4 = g2s(n-4)
          f2snm4 = f2s(n-4)

 4        J = J*q*q
          r = J/a - b

          g1p = g1(n-5) + (3.d0*dx/10.d0)*(11.d0*g1sn - 14.d0*g1snm1 +
     .      26.d0*g1snm2 - 14.d0*g1snm3 + 11.d0*g1snm4)
          f1p = f1(n-5) + (3.d0*dx/10.d0)*(11.d0*f1sn - 14.d0*f1snm1 +
     .      26.d0*f1snm2 - 14.d0*f1snm3 + 11.d0*f1snm4)
          g2p = g2(n-5) + (3.d0*dx/10.d0)*(11.d0*g2sn - 14.d0*g2snm1 +
     .      26.d0*g2snm2 - 14.d0*g2snm3 + 11.d0*g2snm4)
          f2p = f2(n-5) + (3.d0*dx/10.d0)*(11.d0*f2sn - 14.d0*f2snm1 +
     .      26.d0*f2snm2 - 14.d0*f2snm3 + 11.d0*f2snm4)

          nit = 0

 5        g1sng1 = J*(cc*nug1(n+1)*f1p - xk/r*g1p)
          f1snf1 = J*((nuf1(n+1)*g1p - dsqrt(1.d0 - u**2)*bm(n+1)*g2p)/cc
     .  + xk/r*f1p)
          g2sng2 = J*(cc*nug2(n+1)*f2p - xk1/r*g2p)
          f2snf2 = J*((nuf2(n+1)*g2p - dsqrt(1.d0 - u**2)*bm(n+1)*g1p)/cc
     .  + xk1/r*f2p)

          g1c = g1(n-3) + (2.d0*dx/45.d0)*(7.d0*g1sng1 + 32.d0*g1sn +
     .      12.d0*g1snm1 + 32.d0*g1snm2 + 7.d0*g1snm3)
          f1c = f1(n-3) + (2.d0*dx/45.d0)*(7.d0*f1snf1 + 32.d0*f1sn +
     .      12.d0*f1snm1 + 32.d0*f1snm2 + 7.d0*f1snm3)
          g2c = g2(n-3) + (2.d0*dx/45.d0)*(7.d0*g2sng2 + 32.d0*g2sn +
     .      12.d0*g2snm1 + 32.d0*g2snm2 + 7.d0*g2snm3)
          f2c = f2(n-3) + (2.d0*dx/45.d0)*(7.d0*f2snf2 + 32.d0*f2sn +
     .      12.d0*f2snm1 + 32.d0*f2snm2 + 7.d0*f2snm3)

          if(dabs(g1c-g1p) - dabs(g1c)*1.d-12) 36,36,10
 36       if(dabs(f1c-f1p) - dabs(f1c)*1.d-12) 6,6,10
 6        if(dabs(g2c-g2p) - dabs(g2c)*1.d-12) 26,26,10
 26       if(dabs(f2c-f2p) - dabs(f2c)*1.d-12) 11,11,10

C... Check for convergence
 10       if(nit .gt. 100) then
            print*,dabs(g1c-g1p)
            print*,dabs(f1c-f1p)
            print*,dabs(g2c-g2p)
            print*,dabs(f2c-f2p)
            goto 11
          endif

          nit = nit + 1

          g1p = g1c
          f1p = f1c
          g2p = g2c
          f2p = f2c

          goto 5

 11       n = n + 1

          g1(n) = g1c
          f1(n) = f1c
          g1s(n) = g1sng1
          f1s(n) = f1snf1
          g1snm4 = g1snm3
          f1snm4 = f1snm3
          g1snm3 = g1snm2
          f1snm3 = f1snm2
          g1snm2 = g1snm1
          f1snm2 = f1snm1
          g1snm1 = g1sn
          f1snm1 = f1sn
          g1sn = g1sng1
          f1sn = f1snf1


          g2(n) = g2c
          f2(n) = f2c
          g2s(n) = g2sng2
          f2s(n) = f2snf2
          g2snm4 = g2snm3
          f2snm4 = f2snm3
          g2snm3 = g2snm2
          f2snm3 = f2snm2
          g2snm2 = g2snm1
          f2snm2 = f2snm1
          g2snm1 = g2sn
          f2snm1 = f2sn
          g2sn = g2sng2
          f2sn = f2snf2

          if(n-nrmx) 4,12,12

   12     continue

C---Normalization of g and f at the MT sphere------------------
C   and determination of the energy derivatives
C--------------------------------------------------------------

          norm = 0.0d0

          do n = 2, nrmx
            psi(1,n) = g1(n)
            psi(2,n) = f1(n)
            psi(3,n) = g2(n)
            psi(4,n) = f2(n)

            ksi(1,n) = g1(n)
            ksi(2,n) = f1(n)
            ksi(3,n) = g2(n)
            ksi(4,n) = f2(n)
          enddo

          call product2(psi,ksi,a,b,rofi,nr,e1,vm,z,xk,norm)

C----Normalize g, save it for calc. of small par. and orthogonalization-
C       open(78,file='test.dat',status='unknown')
C       print*,'l,mu,alpha=',l,mu,alpha,e,enu
          do n = 2, nrmx
            gn1(alpha,ie,n) = g1(n)/dsqrt(norm)
            gn2(alpha,ie,n) = g2(n)/dsqrt(norm)
            fn1(alpha,ie,n) = f1(n)/dsqrt(norm)
            fn2(alpha,ie,n) = f2(n)/dsqrt(norm)
C         write(78,*) rofi(n),g1(n),f1(n)
          enddo
C       call zprm('test.dat',1,g1,nrmx,nrmx,1)
C       close(78)
C----End of the Energy loop--------------------------------------------
        enddo

C----End of the alpha loop---------------------------------------------
      enddo


C---Orthoganalization--------------------------------------------
      do ie= 1, 5
        e1= enu1 + (ie - 3)*de
        e2= enu2 + (ie - 3)*de
        product = 0.0d0

        do n = 2, nrmx
          psi(1,n) = gn1(1,ie,n)
          psi(2,n) = fn1(1,ie,n)
          psi(3,n) = gn2(1,ie,n)
          psi(4,n) = fn2(1,ie,n)

          ksi(1,n) = gn1(2,ie,n)
          ksi(2,n) = fn1(2,ie,n)
          ksi(3,n) = gn2(2,ie,n)
          ksi(4,n) = fn2(2,ie,n)
C         if (l .eq. 1 .and. mu .eq. -0.5) then
C           write(78,*) rofi(n),gn1(1,3,n),gn2(2,3,n)
C         endif
        enddo

        call product2(psi,ksi,a,b,rofi,nr,enu1,vm,z,xk,product)
cc      write(78,*) 'before', l, product

        do n = 2, nrmx
          gn1(2,ie,n) = -product*gn1(1,ie,n) + gn1(2,ie,n)
          gn2(2,ie,n) = -product*gn2(1,ie,n) + gn2(2,ie,n)
          fn1(2,ie,n) = -product*fn1(1,ie,n) + fn1(2,ie,n)
          fn2(2,ie,n) = -product*fn2(1,ie,n) + fn2(2,ie,n)
        enddo

C---Test orthoganalization

cc      product = 0.0d0

cc      do n = 2, nrmx
cc         ksi(1,n) = gn1(2,ie,n)
cc         ksi(2,n) = fn1(2,ie,n)
cc         ksi(3,n) = gn2(2,ie,n)
cc         ksi(4,n) = fn2(2,ie,n)
cc      enddo

cc      call product2(psi,ksi,a,b,rofi,nr,enu1,vm,z,xk,product)
cc      write(78,*) 'after', l, product


C----Save the ortho-normalized g and f at the MT sphere and calc-------
C----the energy derivative at the enu at the MT sphere          -------

        do alpha = 1, 2
          gmt0(1,alpha,ie) = gn1(alpha,ie,nr)/rofi(nr)
          gmt0(2,alpha,ie) = gn2(alpha,ie,nr)/rofi(nr)
          fmt0(1,alpha,ie) = fn1(alpha,ie,nr)/rofi(nr)
          fmt0(2,alpha,ie) = fn2(alpha,ie,nr)/rofi(nr)
        enddo

C---End of the orthog. energy loop-------------------------------------
      enddo

cc      if(l .eq. 2 ) then
cc         write(78,*) 'mu=', mu
cc         do n = 2, nr
cc         write(78,*) rofi(n), gn1(1,3,n)/rofi(n)
cc         enddo
cc      endif

C-----Save radial derivative at the MT sphere--------------------------


      nug1(nr) = 1.d0 + (1.d0/csq*(enu1 + 2.d0*z/rofi(nr)
     .            - vm(nr) + up*bm(nr)))

C---Added to take care of the abs(mu)=l+0.5 case-----------
      if(dabs(mu) .eq. (l+0.5d0)) then
        nug1(nr) = 1.d0 + (1.d0/csq*(enu2 + 2.d0*z/rofi(nr)
     .            - vm(nr) - upp*bm(nr)))
      endif
C----------------------------------------------------------

      nug2(nr) = 1.d0 + (1.d0/csq*(enu2 + 2.d0*z/rofi(nr)
     .            - vm(nr) - upp*bm(nr)))


      do alpha = 1, 2

        gsmt(1,alpha) = (cc*nug1(nr)*fn1(alpha,3,nr)
     .          - xk/rofi(nr)*gn1(alpha,3,nr))
        gsmt(2,alpha) = (cc*nug2(nr)*fn2(alpha,3,nr)
     .          - xk1/rofi(nr)*gn2(alpha,3,nr))

      enddo


C--Save g and f at the MT sphere--------------------------------------

C ... Multiply the second soltion (alpha=2) by a phase factor (-1)
      if (gmt0(1,1,3) .gt. 0d0 .and. gmt0(2,2,3) .lt. 0d0) then
        print*,'multiplied by a phase factor'
        do i = 1, 2
          gmt0(i,2,3)=-gmt0(i,2,3)
          fmt0(i,2,3)=-fmt0(i,2,3)
        enddo
      endif
      if (gmt0(1,1,3) .lt. 0d0 .and. gmt0(2,2,3) .gt. 0d0) then
        do i = 1, 2
          gmt0(i,2,3)=-gmt0(i,2,3)
          fmt0(i,2,3)=-fmt0(i,2,3)
        enddo
      endif

      do i = 1, 2
        do alpha = 1, 2

          gmt(i,alpha) = gmt0(i,alpha,3)
          fmt(i,alpha) = fmt0(i,alpha,3)

        enddo
      enddo


C-----Calculate the energy derivatives of rg, rf at every point r -----

      do n= 2, nrmx
        do alpha= 1, 2
          gn1dot(alpha,n) = (-8.d0*(gn1(alpha,2,n) - gn1(alpha,4,n))
     .                      + gn1(alpha,1,n) - gn1(alpha,5,n))/12.d0/de
          gn2dot(alpha,n) = (-8.d0*(gn2(alpha,2,n) - gn2(alpha,4,n))
     .                      + gn2(alpha,1,n) - gn2(alpha,5,n))/12.d0/de
          fn1dot(alpha,n) = (-8.d0*(fn1(alpha,2,n) - fn1(alpha,4,n))
     .                      + fn1(alpha,1,n) - fn1(alpha,5,n))/12.d0/de
          fn2dot(alpha,n) = (-8.d0*(fn2(alpha,2,n) - fn2(alpha,4,n))
     .                      + fn2(alpha,1,n) - fn2(alpha,5,n))/12.d0/de
        enddo
      enddo

C--Orthogonality of <phi|phidot>------ same alpha --------------------

      do alpha = 1, 2
        product = 0.0d0
        do n = 2, nrmx
          psi(1,n) = gn1(alpha,3,n)
          psi(2,n) = fn1(alpha,3,n)
          psi(3,n) = gn2(alpha,3,n)
          psi(4,n) = fn2(alpha,3,n)

          ksi(1,n) = gn1dot(alpha,n)
          ksi(2,n) = fn1dot(alpha,n)
          ksi(3,n) = gn2dot(alpha,n)
          ksi(4,n) = fn2dot(alpha,n)
        enddo

        call product2(psi,ksi,a,b,rofi,nr,enu1,vm,z,xk,product)

        do n = 2, nrmx
          gn1dot(alpha,n) = -product*gn1(alpha,3,n) + gn1dot(alpha,n)
          gn2dot(alpha,n) = -product*gn2(alpha,3,n) + gn2dot(alpha,n)
          fn1dot(alpha,n) = -product*fn1(alpha,3,n) + fn1dot(alpha,n)
          fn2dot(alpha,n) = -product*fn2(alpha,3,n) + fn2dot(alpha,n)
        enddo

      enddo

C--Orthogonality of <phi|phidot>------ different alphas --------------

      product = 0.0d0

      do n = 2, nrmx

        psi(1,n) = gn1(1,3,n)
        psi(2,n) = fn1(1,3,n)
        psi(3,n) = gn2(1,3,n)
        psi(4,n) = fn2(1,3,n)

        ksi(1,n) = gn1dot(2,n)
        ksi(2,n) = fn1dot(2,n)
        ksi(3,n) = gn2dot(2,n)
        ksi(4,n) = fn2dot(2,n)

      enddo

      call product2(psi,ksi,a,b,rofi,nr,enu1,vm,z,xk,product)

      product1 = 0.0d0

      do n = 2, nrmx

        psi(1,n) = gn1(2,3,n)
        psi(2,n) = fn1(2,3,n)
        psi(3,n) = gn2(2,3,n)
        psi(4,n) = fn2(2,3,n)

        ksi(1,n) = gn1dot(1,n)
        ksi(2,n) = fn1dot(1,n)
        ksi(3,n) = gn2dot(1,n)
        ksi(4,n) = fn2dot(1,n)

      enddo

      call product2(psi,ksi,a,b,rofi,nr,enu1,vm,z,xk,product1)
cc      write(78,*) 'before', l, product1

      do n = 2, nrmx
        gn1dot(1,n) = -product1*gn1(2,3,n) + gn1dot(1,n)
        gn2dot(1,n) = -product1*gn2(2,3,n) + gn2dot(1,n)
        fn1dot(1,n) = -product1*fn1(2,3,n) + fn1dot(1,n)
        fn2dot(1,n) = -product1*fn2(2,3,n) + fn2dot(1,n)
      enddo

      do n = 2, nrmx
        gn1dot(2,n) = -product*gn1(1,3,n) + gn1dot(2,n)
        gn2dot(2,n) = -product*gn2(1,3,n) + gn2dot(2,n)
        fn1dot(2,n) = -product*fn1(1,3,n) + fn1dot(2,n)
        fn2dot(2,n) = -product*fn2(1,3,n) + fn2dot(2,n)
      enddo

      product = 0.d0
      product1 = 0.d0

      do n = 2, nrmx

        ksi(1,n) = gn1dot(1,n)
        ksi(2,n) = fn1dot(1,n)
        ksi(3,n) = gn2dot(1,n)
        ksi(4,n) = fn2dot(1,n)

      enddo

      call product2(psi,ksi,a,b,rofi,nr,enu1,vm,z,xk,product)
cc       write(78,*) 'after', l, product

C----Save gdot and fdot at the MT sphere-----------------------------

      do alpha = 1, 2
        gmtde(1,alpha) = gn1dot(alpha,nr)/rofi(nr)
        gmtde(2,alpha) = gn2dot(alpha,nr)/rofi(nr)
        fmtde(1,alpha) = fn1dot(alpha,nr)/rofi(nr)
        fmtde(2,alpha) = fn2dot(alpha,nr)/rofi(nr)
      enddo

C-----Calculate the small parameter p=<gdot|gdot>--------------------
C----Try Trapezium rule------------------------------------------------

      do alpha=1, 2
        do i= 1, 2
          smalpa(i,alpha) = 0.0d0
          smalpa1(i,alpha) = 0.0d0
        enddo
      enddo

      enu = 0.5d0*(enum1+enum2)

C      do alpha = 1, 2
C       do n = 2, nrmx-1
C
C       factor(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C       factor1(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n+1)
C     .                - vm(n+1))/cc)*rofi(n+1))**2)
C       factor2(n) = (1d0 + xk1*(xk1+1d0)/((cc + (enu + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C       factor3(n) = (1d0 + xk1*(xk1+1d0)/((cc + (enu + 2d0*z/rofi(n+1)
C     .                - vm(n+1))/cc)*rofi(n+1))**2)



C         smalpa(1,alpha)= smalpa(1,alpha) + 0.5d0*(a*(rofi(n)+b)*
C     .          (gn1dot(alpha,n)*gn1dot(alpha,n)*factor(n) +
C     .           fn1dot(alpha,n)*fn1dot(alpha,n)) +
C     .      a*(rofi(n+1)+b)*(gn1dot(alpha,n+1)*gn1dot(alpha,n+1)*
C     .      factor1(n) + fn1dot(alpha,n+1)*fn1dot(alpha,n+1)))

C         smalpa(2,alpha)= smalpa(2,alpha) +  0.5d0*(a*(rofi(n)+b)*
C     .             (gn2dot(alpha,n)*gn2dot(alpha,n)*factor2(n) +
C     .              fn2dot(alpha,n)*fn2dot(alpha,n)) +
C     .      a*(rofi(n+1)+b)*(gn2dot(alpha,n+1)*gn2dot(alpha,n+1)*
C     .      factor3(n) + fn2dot(alpha,n+1)*fn2dot(alpha,n+1)))

C      enddo
C      enddo

C----Calculate the small parameter p=<Phi^a|Phi^a'>

      do alp1 = 1, 2
        do alp2 = 1, 2

          do n = 2, nrmx-1

            factor(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n)
     .                - vm(n))/cc)*rofi(n))**2)
            factor1(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n+1)
     .                - vm(n+1))/cc)*rofi(n+1))**2)
            factor2(n) = (1d0 + xk1*(xk1+1d0)/((cc + (enu + 2d0*z/rofi(n)
     .                - vm(n))/cc)*rofi(n))**2)
            factor3(n) = (1d0 + xk1*(xk1+1d0)/((cc + (enu + 2d0*z/rofi(n+1)
     .                - vm(n+1))/cc)*rofi(n+1))**2)


            smalpa(alp1,alp2) = smalpa(alp1,alp2) + 0.5d0*(a*(rofi(n)+b)*
     .          (gn1dot(alp1,n)*gn1dot(alp2,n)*factor(n) +
     .           fn1dot(alp1,n)*fn1dot(alp2,n)) +
     .      a*(rofi(n+1)+b)*(gn1dot(alp1,n+1)*gn1dot(alp2,n+1)*
     .      factor1(n) + fn1dot(alp1,n+1)*fn1dot(alp2,n+1)))
     .     +  0.5d0*(a*(rofi(n)+b)*
     .             (gn2dot(alp1,n)*gn2dot(alp2,n)*factor2(n) +
     .              fn2dot(alp1,n)*fn2dot(alp2,n)) +
     .      a*(rofi(n+1)+b)*(gn2dot(alp1,n+1)*gn2dot(alp2,n+1)*
     .      factor3(n) + fn2dot(alp1,n+1)*fn2dot(alp2,n+1)))


          enddo

C       do n=2,nrmx-2
C       factor(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C       factor1(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n+1)
C     .                - vm(n+1))/cc)*rofi(n+1))**2)
C       factor4(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n+2)
C     .                - vm(n+2))/cc)*rofi(n+2))**2)
C      factor2(n) = (1d0 + xk1*(xk1+1d0)/((cc + (enu + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C       factor3(n) = (1d0 + xk1*(xk1+1d0)/((cc + (enu + 2d0*z/rofi(n+1)
C     .                - vm(n+1))/cc)*rofi(n+1))**2)
C       factor5(n) = (1d0 + xk1*(xk1+1d0)/((cc + (enu + 2d0*z/rofi(n+2)
C     .                - vm(n+2))/cc)*rofi(n+2))**2)
C
C       smalpa1(alp1,alp2) = smalpa1(alp1,alp2) + (a*(rofi(n)+b)*
C     .          (gn1dot(alp1,n)*gn1dot(alp2,n)*factor(n) +
C     .           fn1dot(alp1,n)*fn1dot(alp2,n)) + 4.d0*
C     .      a*(rofi(n+1)+b)*(gn1dot(alp1,n+1)*gn1dot(alp2,n+1)*
C     .      factor1(n) + fn1dot(alp1,n+1)*fn1dot(alp2,n+1))+
C     .       a*(rofi(n+2)+b)*(gn1dot(alp1,n+2)*gn1dot(alp2,n+2)*
C     .      factor4(n) + fn1dot(alp1,n+2)*fn1dot(alp2,n+2)))/3.d0
C     .     + (a*(rofi(n)+b)*
C     .             (gn2dot(alp1,n)*gn2dot(alp2,n)*factor2(n) +
C     .              fn2dot(alp1,n)*fn2dot(alp2,n)) + 4.d0*
C     .      a*(rofi(n+1)+b)*(gn2dot(alp1,n+1)*gn2dot(alp2,n+1)*
C     .      factor3(n) + fn2dot(alp1,n+1)*fn2dot(alp2,n+1)) +
C     .         a*(rofi(n+2)+b)*(gn2dot(alp1,n+2)*gn2dot(alp2,n+2)*
C     .      factor5(n) + fn2dot(alp1,n+2)*fn2dot(alp2,n+2)))/3.d0

C          enddo

        enddo
      enddo


C---Try Simpson's Rule-------------------------------------------------

cc      do alpha = 1, 2
cc      do n = 3, nrmx-2, 2

cc          factor(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n)
cc     .                - vm(n))/cc)*rofi(n))**2)

cc         smalpa(1,alpha)= smalpa(1,alpha) + (rofi(n)+b)*
cc     .          (gn1dot(alpha,n)*gn1dot(alpha,n)*factor(n) +
cc     .           fn1dot(alpha,n)*fn1dot(alpha,n))

cc         smalpa(2,alpha)= smalpa(2,alpha) + (rofi(n)+b)*
cc     .             (gn2dot(alpha,n)*gn2dot(alpha,n)*factor(n) +
cc     .           fn2dot(alpha,n)*fn2dot(alpha,n))

cc      enddo
cc      enddo

cc      do alpha = 1, 2
cc       smalpa(1,alpha) = 2.d0*smalpa(1,alpha)
cc       smalpa(2,alpha) = 2.d0*smalpa(2,alpha)
cc      enddo

cc       do alpha = 1, 2
cc       do n = 4, nrmx-1, 2

cc          factor(n) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(n)
cc     .                - vm(n))/cc)*rofi(n))**2)

cc         smalpa(1,alpha)= smalpa(1,alpha) + (rofi(n)+b)*
cc     .          (gn1dot(alpha,n)*gn1dot(alpha,n)*factor(n) +
cc     .           fn1dot(alpha,n)*fn1dot(alpha,n))

cc         smalpa(2,alpha)= smalpa(2,alpha) + (rofi(n)+b)*
cc     .             (gn2dot(alpha,n)*gn2dot(alpha,n)*factor(n) +
cc     .           fn2dot(alpha,n)*fn2dot(alpha,n))

cc      enddo
cc      enddo


cc      factor(2) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(2)
cc     .                - vm(2))/cc)*rofi(2))**2)

cc      factor(nr) = (1d0 + xk*(xk+1d0)/((cc + (enu + 2d0*z/rofi(nr)
cc     .                - vm(nr))/cc)*rofi(nr))**2)

cc      do alpha = 1, 2
cc      smalpa(1,alpha) = (2.d0*smalpa(1,alpha) +
cc     .                   gn1dot(alpha,2)*gn1dot(alpha,2)*factor(2) +
cc     .                   fn1dot(alpha,2)*fn1dot(alpha,2) +
cc     .                   gn1dot(alpha,nr)*gn1dot(alpha,nr)*factor(nr)
cc     .                   + fn1dot(alpha,nr)*fn1dot(alpha,nr))*a/3.d0
cc      smalpa(2,alpha) = (2.d0*smalpa(2,alpha) +
cc     .                   gn2dot(alpha,2)*gn2dot(alpha,2)*factor(2) +
cc     .                   fn2dot(alpha,2)*fn2dot(alpha,2) +
cc     .                   gn2dot(alpha,nr)*gn2dot(alpha,nr)*factor(nr)
cc     .                   + fn2dot(alpha,nr)*fn2dot(alpha,nr))*a/3.d0

cc      enddo

cc      if(dabs(mu) .eq. (l+0.5d0)) then
cc       write(78,*) mu, smalpa(1,1)
cc      endif

cc       enddo
cc      enddo


C ... Copy small parameter to pprel(4)
C     open(180,file='pa.dat',status='unknown')
      do i = 1, 2
        do k = 1, 2
          pprel(4,l,imu,i,k) = smalpa(i,k)
C         write(180,*) smalpa(i,k)
        enddo
      enddo

Cgetarg       end
      end subroutine rdeq


      subroutine fdpp(enu1,enu2,ksop,shft,gmt,fmt,gmtde,fmtde,z,rmax,
     .  avw,l,lmx,imu,sr1,sr2,srav1,srav2,srdot1,srdot2,sravdot1,
     .  sravdot2,gsr1,gsr2,gsrav1,gsrav2,gsrdot1,gsrdot2,gsravdot1,
     .  gsravdot2,pprel,gsmt)

C------------------------------------------------------------------------------
C      Makes potential parameters
C      Inputs:
C      gmt,fmt: 2x2 normalized wave function x r, at the MT surface
C      gmtde, fmtde: d/de(gmt) and d/de(fmt)
C      l, mu: quantum numbers
C-------------------------------------------------------------------------------
C     implicit none
C     Input parameters
      integer l,lmx,imu
      double precision gmt(2,2), fmt(2,2), gmtde(2,2), fmtde(2,2),
     .                 mu, rmax, avw, z, enu1, enu2, sr1, sr2, srav1,
     .                 srav2, srdot1, srdot2, sravdot1, sravdot2,
     .                 gsrav1, gsrav2, gsrdot1, gsrdot2, gsravdot1,
     .                 gsravdot2, gsr1, gsr2,gsmt(2,2),shft
C     Local variables
      integer  i, j
C     double precision ftmt(2,2),fmtdet(2,2),temp(2,2)
      double precision cc, gtmt(2,2), gmt1(2,2),
     .                 gmtde1(2,2),fg(2,2),
     .                 fgde(2,2), xk, xk1, kapa(2,2), D(2,2),
     .                 Ddot(2,2), Dtdot(2,2), lmatrix(2,2),
     .                 delta(2,2), delta1(2,2), deltat(2,2),
     .                 deltat1(2,2), unit(2,2), q(2,2), q1(2,2),
     .                 gam(2,2), gam1(2,2), aa(2,2), aa1(2,2),
     .                 one(2,2), two(2,2), gmtdet(2,2), gamt(2,2),
     .                 yv(2,2), yv1(2,2), vm(2,2), enum(2,2),
     .                 gamma(2,2), gammat(2,2), delt(2,2),
     .                 deltt(2,2), cm(2,2), cm1(2,2),
     .                 fg1(2,2), savw, ksop(0:lmx,2,2,6)
C     double precision clebsh1(2,2), clebsh1t(2,2),clebsh(2,2), u
C     Output
      double precision pprel(4,0:lmx,2*(lmx+1),2,2)
c     double precision pprel(4,2*(lmx+1),2,2)

      if(z .eq. 0) goto 210

      mu = dble(imu-l) - 1.5d0
      cc = 274.072d0
      savw = (avw/rmax)**(2*l+1)
cc------Test:Eliminate Spin-Orbit for non-mag. materials
cc      mu = dble(l)+0.5d0

C-------Definition of useful matrices for the subroutine-----

      if (dabs(mu) .eq. (l+0.5d0)) then
        xk1= dble(-l)
        xk = xk1
      else
        xk1=dble(-l)
        xk =dble(l+1)
      endif

      unit(1,1) = 1.d0
      unit(2,2) = 1.d0
      unit(1,2) = 0.0d0
      unit(2,1) = 0.0d0

c         print*,'-------------------------'
c         print*, l,ksop(l,1,1,1),enu1,enu2

      enum(1,1) = -0.5d0*dble((l+1))*ksop(l,1,1,1)
     .               + 0.5d0*(enu1 + enu2)
      enum(2,2) = 0.5d0*dble(l)*ksop(l,1,1,1)
     .               + 0.5d0*(enu1 + enu2)

      if (dabs(mu) .eq. (l+0.5d0)) enum(1,1) = enum(2,2)

      enum(1,2) = 0.0d0
      enum(2,1) = 0.0d0

c         print*, l,enum(1,1),enum(2,2)


      kapa(1,1) = xk
      kapa(2,2) = xk1
      kapa(1,2) = 0.0d0
      kapa(2,1) = 0.0d0

      lmatrix(1,1) = dble(l)
      lmatrix(2,2) = dble(l)
      lmatrix(1,2) = 0.0d0
      lmatrix(2,1) = 0.0d0

cc        if(l .eq. 2) then
cc         print*, '===================='
cc         print*, gmt(1,1), gmt(1,2)
cc         print*, gmt(2,1), gmt(2,2)
cc        endif

C------Invert g and g_dot----------------------------
C      The g and f do not commute. In Phys Rev. B
C      43, 14414, 1991 it is D = c*S*g^(-1)*f - k - I
C      In Sov. Phys. Solid State 31(8), 1285 1989 is
C      D = c*S*f*g^(-1) - k - I. However only the second
C      formula gives symmetric D. The order will depend
C      on the order in g and f. If alpha numerates columns then
C      the second formula works, if alpha numerates strokes then
C      the first formula should be applied.
C----------------------------------------------------

      call dinv22(gmt,gmt1)
      call dinv22(gmtde,gmtde1)

      call mul22(fmt,gmt1,fg)
      call mul22(fmtde,gmtde1,fgde)

      call mul22(gsmt,gmt1,fg1)

C----Logarithmic derivative matrices-----------------

      do i = 1, 2
        do j = 1, 2
          D(i,j) = rmax*cc*fg(i,j) - kapa(i,j)
c           D(i,j) = fg1(i,j) - unit(i,j)
          Ddot(i,j) = rmax*cc*fgde(i,j) - kapa(i,j)
        enddo
      enddo

C------Symmetrize D (it is symmetric, but for numerical)

      D(1,2) = (D(1,2) + D(2,1))/2.d0
      D(2,1) = D(1,2)

C------Make corrections for diagonal D, D_dot--------

cc       if(mu .lt. 0) then
cc       D(1,1) = D(1,1) - srav1 + sr1
cc       D(2,2) = D(2,2) - srav1 + sr1
cc       Ddot(1,1) = Ddot(1,1) - sravdot1 + srdot1
cc       Ddot(2,2) = Ddot(2,2) - sravdot1 + srdot1
cc       else
cc       D(1,1) = D(1,1) - srav2 + sr2
cc       D(2,2) = D(2,2) - srav2 + sr2
cc       Ddot(1,1) = Ddot(1,1) - sravdot2 + srdot2
cc       Ddot(2,2) = Ddot(2,2) - sravdot2 + srdot2
cc       endif

C------Make corrections for diagonal g, g_dot--------

cc       if(mu .lt. 0) then
cc       gmt(1,1) = gmt(1,1) - gsrav1 + gsr1
cc       gmt(2,2) = gmt(2,2) - gsrav1 + gsr1
cc       gmtde(1,1) = gmtde(1,1) - gsravdot1 + gsrdot1
cc       gmtde(2,2) = gmtde(2,2) - gsravdot1 + gsrdot1
cc       else
cc       gmt(1,1) = gmt(1,1) - gsrav2 + gsr2
cc       gmt(2,2) = gmt(2,2) - gsrav2 + gsr2
cc       gmtde(1,1) = gmtde(1,1) - gsravdot2 + gsrdot2
cc       gmtde(2,2) = gmtde(2,2) - gsravdot2 + gsrdot2
cc      endif

C--Transpose gmt and Ddot -------------------------

      do i = 1, 2
        do j = 1, 2
C            ftmt(i,j) = fmt(j,i)
C            fmtdet(i,j) = fmtde(j,i)
          gtmt(i,j) = gmt(j,i)
          gmtdet(i,j) = gmtde(j,i)
          Dtdot(i,j) = Ddot(j,i)
        enddo
      enddo

C---Make unscreened potential parameters-------------

      do i = 1, 2
        do j = 1, 2
          delta(i,j) = D(i,j) - Ddot(i,j)
          deltat(i,j) = D(i,j) - Dtdot(i,j)
        enddo
      enddo

      call dinv22(delta, delta1)
      call dinv22(deltat, deltat1)

cc       call mul22(gmtdet,gmt,one)
cc       call mul22(gtmt,gmtde,two)

      call mul22(gmt,gmtdet,one)
      call mul22(gmtde,gtmt,two)

      do i = 1, 2
        do j = 1, 2
          aa(i,j) = -0.5d0*(delta1(i,j)+deltat1(i,j))
cc        aa(i,j) = -0.5d0*rmax*(one(i,j) + two(i,j))
        enddo
      enddo

      call dinv22(aa,aa1)

      do i = 1, 2
        do j = 1, 2
          delta(i,j) = D(i,j) - lmatrix(i,j)
        enddo
      enddo

      call dinv22(delta,delta1)

      do i = 1, 2
        do j = 1, 2
          q(i,j) = delta(i,j) + aa1(i,j)
        enddo
      enddo

      call dinv22(q,q1)

      do i = 1, 2
        do j = 1, 2
          q(i,j) = 2*(2*l+1)*savw*((2*l+1)*q1(i,j)+unit(i,j))
        enddo
      enddo

      call mul22(aa,delta,gam)

      do i = 1, 2
        do j = 1, 2
          gam(i,j) = gam(i,j) + unit(i,j)
        enddo
      enddo

      call dinv22(gam,gam1)
cc       call mul22(gtmt,gam1,gam)
      call mul22(gam1,gmt,gam)


      do i = 1, 2
        do j = 1, 2
          gam(i,j) = dsqrt(2.d0*rmax)*dsqrt(savw)*(2*l+1)*gam(i,j)
        enddo
      enddo

      do i = 1, 2
        do j = 1, 2
          gamt(i,j) = gam(j,i)
        enddo
      enddo

      do i = 1, 2
        do j = 1, 2
          yv(i,j) = delta1(i,j) + aa(i,j)
        enddo
      enddo

      call dinv22(yv,yv1)
      call mul22(gtmt,yv1,yv)
      call mul22(yv,gmt,vm)

      do i = 1, 2
        do j = 1, 2
          vm(i,j) = rmax*vm(i,j) + enum(i,j)
        enddo
      enddo

C-----------TEST----------------------------------------------------

cc       do i = 1, 2
cc          do j = 1, 2
cc            delta(i,j) = D(i,j) + lmatrix(i,j) + unit(i,j)
cc           enddo
cc       enddo

cc       call dinv22(delta,delta1)

cc        do i = 1, 2
cc          do j = 1, 2
cc         yv(i,j) = delta1(i,j) + aa(i,j)
cc          enddo
cc       enddo

cc       call dinv22(yv,yv1)
cc       call mul22(gtmt,yv1,temp)
cc       call mul22(temp,gmt,cm)

cc       do i = 1, 2
cc          do j = 1, 2
cc             cm(i,j) = rmax*cm(i,j) + enum(i,j)
cc         enddo
cc      enddo

cc       do i = 1, 2
cc          do j = 1, 2
cc             delta(i,j) = Ddot(i,j) + lmatrix(i,j) + unit(i,j)
cc             delta1(i,j) = Ddot(i,j) - lmatrix(i,j)
cc          enddo
cc       enddo

cc       call dinv22(delta1,temp)
cc       call mul22(delta,temp,q1)

cc       do i = 1, 2
cc          do j = 1, 2
cc             q(i,j) = 2*(2*l+1)*savw*q1(i,j)
cc          enddo
cc       enddo


C-----Potential parameters in the nearly orthonormalized-----
C     representation
C------------------------------------------------------------

      call dinv22(q,gamma)

      do i = 1, 2
        do j = 1, 2
          gammat(i,j) = gamma(j,i)
        enddo
      enddo

cc       call mul22(gam,gamma,delt)
cc       call mul22(gammat,gamt,deltt)

      call mul22(gamma,gam,delt)
      call mul22(gamt,gamma,deltt)

cc      call mul22(gam,gamma,cm)
cc     call mul22(cm,gamt,cm1)

cc       call mul22(delt,q,cm)
cc       call mul22(cm,deltt,cm1)

      call mul22(deltt,q,cm)
      call mul22(cm,delt,cm1)


      do i = 1, 2
        do j = 1, 2
          cm(i,j) = vm(i,j) + cm1(i,j)
        enddo
        cm(i,i) = cm(i,i) + shft
      enddo

C       open(190,file='cm.dat',status='unknown')
C       open(191,file='gamma.dat',status='unknown')
C       open(192,file='delt.dat',status='unknown')

      do i = 1, 2
        do j = 1, 2
          pprel(1,l,imu,i,j) = cm(i,j)
C          write(190,*) cm(i,j)
          pprel(2,l,imu,i,j) = gamma(i,j)
C          write(191,*) gamma(i,j)
          pprel(3,l,imu,i,j) = delt(i,j)
C          write(192,*) delt(i,j)
        enddo
      enddo

C------rotate cm---------------------------------------

cc       u = mu/(dble(l)+0.5d0)

cc       clebsh(1,1) = dsqrt(1.d0+u)/dsqrt(2.d0)
cc       clebsh(2,2) = clebsh(1,1)
cc       clebsh(2,1) = dsqrt(1.d0-u)/dsqrt(2.d0)
cc       clebsh(1,2) = -clebsh(2,1)

cc       call dinv22(clebsh,clebsh1)

cc      call mul22(clebsh1t,cm,cm1)
cc      call mul22(cm1,clebsh1,cm)

cc      if(dabs(mu) .eq. (l+0.5d0)) then
cc      write(77,*) mu, xk-1, xk1-1
cc      write(77,*) '--------------------------------------------'
cc      write(77,*) cm(1,1), cm(1,2)
cc      write(77,*) cm(2,1), cm(2,2)
cc      write(77,*) '--------------------------------------------'
cc      endif


C-------This is a failed way to built the pot. pars-----------
cc       do i = 1, 2
cc          do j = 1, 2
cc             delta(i,j) = D(i,j) - Ddot(i,j)
cc          enddo
cc       enddo

cc       do i = 1, 2
cc          do j = 1, 2
cc             deltat(i,j) = delta(j,i)
cc          enddo
cc       enddo


cc       call dinv22(delta, delta1)
cc       call dinv22(deltat, deltat1)

cc       call mul22(gmtdet,gmt,one)
cc       call mul22(gtmt,gmtde,two)

cc       do i = 1, 2
cc        do j = 1, 2
cc        aa(i,j) = -0.5d0*(delta1(i,j) + deltat1(i,j))
cc         aa(i,j) = -0.5d0*rmax*(one(i,j) + two(i,j))
cc        enddo
cc       enddo

cc       call dinv22(aa,aa1)

cc       do i = 1, 2
cc          do j = 1, 2
cc             delta(i,j) = D(i,j) - lmatrix(i,j)
cc           enddo
cc      enddo

cc       call dinv22(delta,delta1)

cc       do i = 1, 2
cc          do j = 1, 2
cc             q(i,j) = delta(i,j) + aa1(i,j)
cc          enddo
cc       enddo

cc      write(77,*) '--------------------------------------------'
cc      write(77,*) z, l, mu, xk, xk1
cc      write(77,*) '--------------------------------------------'
cc      write(77,*) delta(1,1), delta(1,2)
cc      write(77,*) delta(2,1), delta(2,2)
cc      write(77,*) '--------------------------------------------'
cc      write(77,*)  aa1(1,1), aa1(1,2)
cc      write(77,*)  aa1(2,1), aa1(2,2)

cc       call dinv22(q,q1)

cc       do i = 1, 2
cc          do j = 1, 2
cc             q(i,j) = 2*(2*l+1)*savw*((2*l+1)*q1(i,j)+unit(i,j))
cc          enddo
cc       enddo

cc       if(dabs(mu) .eq. (l+0.5d0)) then
cc          write(77,*) l,q(1,1)
cc       endif

cc       call mul22(aa,delta,gam)

cc       do i = 1, 2
cc          do j = 1, 2
cc             gam(i,j) = gam(i,j) + unit(i,j)
cc          enddo
cc       enddo

cc       call dinv22(gam,gam1)
cc       call mul22(gam1,gmt,gam)

cc       do i = 1, 2
cc          do j = 1, 2
cc             gam(i,j) = dsqrt(2.d0*rmax*savw)*(2*l+1)*gam(i,j)
cc          enddo
cc       enddo

cc       do i = 1, 2
cc          do j = 1, 2
cc             gamt(i,j) = gam(j,i)
cc          enddo
cc       enddo

cc       do i = 1, 2
cc          do j = 1, 2
cc          yv(i,j) = delta1(i,j) + aa(i,j)
cc          enddo
cc       enddo

cc       call dinv22(yv,yv1)
cc       call mul22(gtmt,yv1,yv)
cc       call mul22(yv,gmt,vm)

cc       do i = 1, 2
cc          do j = 1, 2
cc            vm(i,j) = rmax*vm(i,j) + enum(i,j)
cc          enddo
cc       enddo

C-----Potential parameters in the nearly orthonormalized-----
C     representation
C------------------------------------------------------------

cc       call dinv22(q,gamma)

cc       do i = 1, 2
cc          do j = 1, 2
cc             gammat(i,j) = gamma(j,i)
cc          enddo
cc       enddo

cc       call mul22(gam,gamma,delt)

cc       call mul22(gam,gamma,cm)
cc       call mul22(cm,gamt,cm1)

cc       do i = 1, 2
cc          do j = 1, 2
cc            cm(i,j) = vm(i,j) + cm1(i,j)
cc          enddo
cc       enddo

cc       open(190,file='cm.dat',status='unknown')
cc       open(191,file='gamma.dat',status='unknown')
cc       open(192,file='delt.dat',status='unknown')

cc       do i = 1, 2
cc        do j = 1, 2
cc           pprel(1,l,imu,i,j) = cm(i,j)
cc           write(190,*) cm(i,j)
cc           pprel(2,l,imu,i,j) = gamma(i,j)
cc           write(191,*) gamma(i,j)
cc           pprel(3,l,imu,i,j) = delt(i,j)
cc           write(192,*) delt(i,j)
cc        enddo
cc       enddo

C-------End of the failed way--------------------------------------

 210  return

Cgetarg        end
      end subroutine fdpp


      subroutine mul22(a,b,c)

C      Multiplies 2x2 matrices a*b
C      a is from the left be carefull for non commuting matrices
C     implicit none
      double precision a(2,2), b(2,2), c(2,2)

      c(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
      c(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2)
      c(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
      c(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2)

      return
Cgetarg        end
      end subroutine mul22



      subroutine product2(psi,ksi,a,b,rofi,nr,e,vm,z,xk,product)

C     Finds the product of two 4-vectors <psi|ksi> using the
C     Trapezium rule
C     Input parameters
      integer nr
      double precision psi(4,nr), ksi(4,nr), a, b, xk, rofi(nr), e,
     .                 vm(nr), z
C     Output
      double precision product
C     Local
      integer n, m
C     double precision cc

C     cc = 274.072d0
      product = 0.0d0

      do n = 2, nr-1
        m = n + 1

C         factor = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C         factor1 = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(m)
C     .                - vm(m))/cc)*rofi(m))**2)

        product = product +  0.5d0*a*((rofi(n)+b)*
     .          ((psi(1,n)*ksi(1,n) + psi(3,n)*ksi(3,n)) +
     .            psi(2,n)*ksi(2,n) + psi(4,n)*ksi(4,n)) +
     .         (rofi(m)+b)*
     .          ((psi(1,m)*ksi(1,m) + psi(3,m)*ksi(3,m)) +
     .            psi(2,m)*ksi(2,m) + psi(4,m)*ksi(4,m)))

      enddo
Cgetarg       end
      end subroutine product2


C      subroutine product4(psi,ksi,a,b,rofi,nr,e,vm,z,xk,product)
CC     Finds the product of two 4-vectors <psi|ksi> using the
CC     Simpson's rule
CC     Input parameters
C      integer nr
C      double precision psi(4,nr), ksi(4,nr), a, b, xk, rofi(nr), e,
C     .                 vm(nr), z
CC     Output
C      double precision product, product1, product2
CC     Local
C      integer n, m
C      double precision factor, factor1, cc
C
C      cc = 274.072d0
C      product = 0.0d0
C      product2 = 0.0d0
C      product1 = 0.0d0
C
C
C      do n = 3, nr-2, 2
C
C          factor = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C
C         product = product + (rofi(n)+b)*
C     .            ((psi(1,n)*ksi(1,n) + psi(3,n)*ksi(3,n))*factor +
C     .            psi(2,n)*ksi(2,n) + psi(4,n)*ksi(4,n))
C
C
C      enddo
C
C      product = (2.0d0*product)
C
C      do n = 4, nr-1, 2
C
C         factor = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C
C         product = product + (rofi(n)+b)*
C     .            ((psi(1,n)*ksi(1,n) + psi(3,n)*ksi(3,n))*factor +
C     .            psi(2,n)*ksi(2,n) + psi(4,n)*ksi(4,n))
C
C         enddo
C
C
C      factor = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(2)
C     .                - vm(2))/cc)*rofi(2))**2)
C
C      factor1 = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(nr)
C     .                - vm(nr))/cc)*rofi(nr))**2)
C
C      product1 = (rofi(2)+b)*
C     .            ((psi(1,2)*ksi(1,2) + psi(3,2)*ksi(3,2))*factor +
C     .            psi(2,n)*ksi(2,n) + psi(4,2)*ksi(4,2))
C
C      product2 = (rofi(nr)+b)*
C     .            ((psi(1,nr)*ksi(1,nr) + psi(3,nr)*ksi(3,nr))*factor1 +
C     .            psi(2,nr)*ksi(2,nr) + psi(4,nr)*ksi(4,nr))
C
C
C      product = (2.0d0*product + product1 + product2)*a/3.0d0
C
C      return
C      end

      subroutine clebsh(l,imu,a,a1)

C     implicit none
      double precision a(2,2),a1(2,2),cl(2,2),mu,u
C     double precision cl1(2,2)
      integer ms1,ms2,l,imu


      mu= dble(imu - l) - 1.5d0

      u = mu/(dble(l)+0.5d0)

      cl(1,1) = dsqrt(0.5d0*(1.d0+u))
      cl(2,2) = cl(1,1)
      cl(2,1) = -dsqrt(0.5d0*(1.d0-u))
      cl(1,2) = -cl(2,1)

C         call dinv22(cl,cl1)

      do ms1 = 1, 2
        do ms2 = 1, 2

          a1(ms1,ms2) =
     .            cl(1,ms1)*a(1,1)*cl(1,ms2) +
     .            cl(1,ms1)*a(1,2)*cl(2,ms2) +
     .            cl(2,ms1)*a(2,1)*cl(1,ms2) +
     .            cl(2,ms1)*a(2,2)*cl(2,ms2)
        enddo
      enddo

C       --- Test rotate back to lms rep'sn from kappa-mu rep'sn  ---
C          do ms1 = 1, 2
C             do ms2 = 1, 2
C                a(ms1,ms2) =
C     .            cl1(1,ms1)*a1(1,1)*cl1(1,ms2) +
C     .            cl1(1,ms1)*a1(1,2)*cl1(2,ms2) +
C     .            cl1(2,ms1)*a1(2,1)*cl1(1,ms2) +
C     .            cl1(2,ms1)*a1(2,2)*cl1(2,ms2)
C
C              else
C                reslms(m1,m2,ms1,ms2) = 0d0
C              endif
C
C            enddo
C            enddo

      return

Cgetarg           end
      end subroutine clebsh














