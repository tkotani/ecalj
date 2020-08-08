      subroutine a2rotm(strn,linv,ipr,rotm)
C- Generate Euler angles from a sequence of rotations, ascii input
Ci strn: Sequence of rotation matrices.  Syntax: rot1[,rot2...],
Ci   where each rotj looks like (x,y,z)angle.  x:, y: or z: may be
Ci   substituted for (x,y,z).  Example:  strn that makes rotation matrix
Ci   corresponding to the Euler angles alpha=pi/4 beta=pi/3 gamma=pi/2:
Ci   (0,0,1)pi/4,(0,1,0)pi/3,(0,0,1)pi/2   or just  z:pi/4,y:pi/3,z:pi/2
Ci linv: T, return transpose of rotm, corresponding to rotation of
Ci   positions rather than the coordinates (default).
Ci ipr: verbosity
Co rotm: rotation matrix.  Follow with a call to rm2eua to find Euler
Co   angles.
C     implicit none
      logical linv
      character*(*) strn, rchr*7, ss*80
      double precision rotm(3,3)
      double precision alpha,beta,gamma,rotj(3,3),r2(3,3),angle,
     .r(3),pi,ddot,fuzz,ct,theta,phi

      integer ls,itrm,ich,a2vec,ix(3),i,j,irot,ipr,i1mach
      parameter (fuzz = 1d-12)
      data rchr /'(XxYyZz'/


      pi = 4*datan(1d0)
      call dpzero(rotm,9)
      rotm(1,1) = 1
      rotm(2,2) = 1
      rotm(3,3) = 1
      ss = strn
      ls = len(ss)
      irot = 0
      ich = 0
C --- Entry point to accumulate new rotation ---
    5 continue
      irot = irot+1
      call chrps2(ss,rchr,len(rchr),ich,ich,itrm)
      if (itrm .eq. 0) goto 999
      if (itrm .gt. 1) then
        call dpzero(r,3)
        if (ss(ich+2:ich+2) .ne. ':') goto 999
        itrm = itrm/2
        if (itrm .gt. 3) call rx('bug in a2rotm')
   11   r(itrm) = 1
        ich = ich+2
      else
        ich = ich+1
        if (a2vec(ss,ls,ich,4,',)',2,2,3,ix,r) .ne. 3) goto 999
      endif
      if (a2vec(ss,ls,ich,4,', ',2,2,1,ix,angle) .ne. 1) goto 999
      call dpzero(rotj,9)
C ... Rotation about r
      call dscal(3,1/dsqrt(ddot(3,r,1,r,1)),r,1)
      ct = r(3)
      phi = 0
      if ( ct .gt.  1 .and.  ct .lt. 1-fuzz) ct = 1
      if (-ct .gt.  1 .and. -ct .lt. 1-fuzz) ct = -1
      if (dabs(ct) .lt. 1) phi = datan2(r(2),r(1))
      theta = dacos(ct)
      if (theta .lt. 0) theta = 2*pi - theta
      call rotma(phi,theta,angle,rotj)
C ... rotm <- rotj*rotm
      call dpcopy(rotm,r2,1,9,1d0)
      call dmpy(rotj,3,1,r2,3,1,rotm,3,1,3,3,3)
      if (ipr .gt. 40) then
        call rm2eua(rotm,alpha,beta,gamma)
        call awrit4(' a2rotm:  rotation %i, alpha = %1;6d  '//
     .  'beta = %1;6d  gamma = %1;6d',' ',80,i1mach(2),
     .  irot,alpha,beta,gamma)
        print 335, ((rotm(i,j),j=1,3),i=1,3)
  335   format((3f15.9))

C ... code to make beta, phi om from Euler angles
C        st = dsqrt(1-ct**2)
C        print *, 'beta=',2d0*dasin(st*dsin(angle))
C        print *, 'theta,phi=',theta,phi
C        print *, 'phi=',phi,(alpha-gamma-pi)/2
C        if (dabs((phi-(alpha-gamma-pi)/2)) .gt. 1d-10) then
C          print *, 'warning: phi off', need fix this
C          phi = phi - pi
C        endif
C        if (dabs((phi-(alpha-gamma-pi)/2)) .gt. 1d-10)  stop 'phi off'
C        gp = (alpha+gamma)/2
C        print *, 'gp=',gp
C        print *, 'angle=', angle,
C     .    dasin(dsqrt((datan(gp)**2+dsin(beta/2)**2)/(datan(gp)**2+1)))
C        tg = dcos(theta)*dtan(angle)
C        tg = dtan(gp)
C        print *, dcos(angle)**2*tg**2+dsin(beta/2)**2-dsin(angle)**2,
C     .    dtan(gp),tg
C        print *, tg**2+dsin(beta/2)**2-dsin(angle)**2*(1+tg**2)
C        print *, (tg**2+dsin(beta/2)**2)/(1+tg**2)-dsin(angle)**2
C        print *, dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2))-dsin(angle)
C        print *, 'angle=', angle,
C     .    dasin(dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2)))
CC ...   sin and cos angle, each undetermined up to a sign
C        sw = dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2))
C        cw = dsqrt(1-sw**2)
C        print *, 'sin angle=', sin(angle),sw
C        if (tg .lt. 0) cw = -cw
CC ...   theta, assuming ??? cw positive
C        th = datan2(dsin(beta/2),cw*tg)
C        print *, 'theta=', theta, datan2(dsin(beta/2),cw*tg)
C        print *, 'theta=', theta, datan2(dsin(beta/2),-cw*tg)
C        st = dsin(th)
C        ct = dcos(th)
C        om = datan2(tg,ct)
C
Cc ...   check
C        print *, 'ck tg', tg,tan(om)*ct
C        print *, 'ck sb', sin(beta/2),sin(th)*sin(om)
C
CC ok for ROT=(.1,.2,.3).9*pi
CC but backwards for ROT=(.1,.2,-.3)-.9*pi ROT=(.1,-.2,-.3)-.9*pi
C        stop
C
C        call euler(r(1),r(2),r(3),angle,alpha,beta,gamma)

      endif

C --- Pick up another rotation matrix, or quit if done ---
      if (ix(1) .eq. 1) goto 5
      if (linv) then
        call dpcopy(rotm,r2,1,9,1d0)
        do  20  i = 1, 3
          do  20  j = 1, 3
   20   rotm(i,j) = r2(j,i)
        if (ipr .ge. 50) then
          call rm2eua(rotm,alpha,beta,gamma)
          call awrit3(' a2rotm:  inverse, alpha = %1;6d  '//
     .    'beta = %1;6d  gamma = %1;6d',' ',80,i1mach(2),
     .    alpha,beta,gamma)
          print 335, ((rotm(i,j),j=1,3),i=1,3)
        endif
      endif
      return

  999 call rx('a2rotm: rotation must be one of x:#,y:#,z:#,(vec)#')
      end
      subroutine eua2rm(alpha,beta,gamma,r)
C- Generate the rotation matrix corresponding to Euler angles
C ----------------------------------------------------------------------
C  Definition: R = Rz(gamma)*Ry(beta)*Rz(alpha)
Cr This definition has the property that a vector v rotated, i.e.
Cr    R v = v'
Cr is rotated to v' = zhat 
Cr when v is defined by the polar coordinates (alpha,beta)
Cr    vx = cos(alpha) sin(beta)
Cr    vy = sin(alpha) sin(beta)
Cr    vx = cos(beta)
Cr Note also that NB: r^-1 = r+
C ----------------------------------------------------------------------
C     implicit none
      double precision r(3,3),alpha,beta,gamma
C     integer i,j
      double precision ca,cb,cg,sa,sb,sg

      ca = dcos(alpha)
      sa = dsin(alpha)
      cb = dcos(beta)
      sb = dsin(beta)
      cg = dcos(gamma)
      sg = dsin(gamma)

C     Euler angles: R3(=R_z(gamma)) R2(=R_y(beta)) R1(=R_z(alpha))
      r(1,1) =  cg*cb*ca - sg*sa
      r(2,1) = -sg*cb*ca - cg*sa
      r(3,1) =  sb*ca
      r(1,2) =  cg*cb*sa + sg*ca
      r(2,2) = -sg*cb*sa + cg*ca
      r(3,2) =  sb*sa
      r(1,3) = -cg*sb
      r(2,3) =  sg*sb
      r(3,3) =  cb

C      print 335, ((r(i,j),j=1,3),i=1,3)
C  335 format((3f15.9))
      end
      subroutine rm2eua(r,alpha,beta,gamma)
C- Generate Euler angles from a rotation matrix
C     implicit none
      double precision r(3,3),alpha,beta,gamma,fuzz,pi
      integer i1mach,i,j
      double precision ca,cb,cg,sa,sb,sg,amg,apg,xx,yy,d1mach
      parameter (fuzz = 1d-8)

      pi = 4*datan(1d0)
      cb = r(3,3)
      if (dabs(cb) .gt. 1 .and.  dabs(cb) .lt. 1+fuzz) cb = sign(1d0,cb)
C --- Handle case beta is near zero or pi ---
      if (dabs(dabs(cb)-1) .lt. dsqrt(fuzz)) then
        xx = r(1,3)**2 + r(2,3)**2
        yy = r(3,1)**2 + r(3,2)**2
        beta = datan2(dsqrt((xx+yy)/2),r(3,3))
      else
        beta = dacos(cb)
      endif
      if (dsqrt(r(3,2)**2 + r(3,1)**2) .lt. fuzz .or.
     .dsqrt(r(2,3)**2 + r(1,3)**2) .lt. fuzz) then
        if (cb .gt. 0) then
          apg = datan2(r(1,2),r(2,2))
          if (r(3,1)**2+r(3,2)**2 .ne. 0d0) then
            alpha = datan2(r(3,2),r(3,1))
            gamma = apg - alpha
          elseif (r(1,3)**2+r(2,3)**2 .ne. 0d0) then
            gamma = datan2(r(2,3),-r(1,3))
            alpha = apg - gamma
          else
            gamma = 0d0
            alpha = apg
          endif
        else
          amg = datan2(-r(1,2),r(2,2))
          if (r(3,1)**2+r(3,2)**2 .ne. 0d0) then
            alpha = datan2(r(3,2),r(3,1))
            if (dsqrt(r(3,1)**2+r(3,2)**2) .lt. d1mach(3)*2) alpha = 0
            gamma = alpha - amg
          elseif (r(1,3)**2+r(2,3)**2 .ne. 0d0) then
            gamma = datan2(r(2,3),-r(1,3))
            alpha = amg + gamma
          else
            gamma = 0
            alpha = amg
          endif
        endif
      else
        alpha = datan2(r(3,2),r(3,1))
        gamma = datan2(r(2,3),-r(1,3))
      endif
      if (alpha .gt.  pi) alpha = alpha - 2*pi
      if (alpha .lt. -pi) alpha = alpha + 2*pi
      if (gamma .gt.  pi) gamma = gamma - 2*pi
      if (gamma .lt. -pi) gamma = gamma + 2*pi
C      call awrit3(' rm2eua:  alpha = %1;6d  beta = %1;6d'//
C     .  '  gamma = %1;6d',' ',80,i1mach(2),alpha,beta,gamma)
C ... Check validity of entire rotation matrix
      ca = dcos(alpha)
      sa = dsin(alpha)
      cb = dcos(beta)
      sb = dsin(beta)
      cg = dcos(gamma)
      sg = dsin(gamma)
C      print *,   dabs(r(1,1)-(ca*cb*cg-sa*sg)) .gt. fuzz*10 ,
C     .    dabs(r(1,2)-(sa*cb*cg+ca*sg)) .gt. fuzz*10 ,
C     .    dabs(r(1,3)-(-sb*cg)) .gt. fuzz*10 ,
C     .    dabs(r(2,1)-(-ca*cb*sg-sa*cg)) .gt. fuzz*10 ,
C     .    dabs(r(2,2)-(-sa*cb*sg+ca*cg)) .gt. fuzz*10 ,
C     .    dabs(r(2,3)-(sb*sg)) .gt. fuzz*10 ,
C     .    dabs(r(3,1)-(ca*sb)) .gt. fuzz*10 ,
C     .    dabs(r(3,2)-(sa*sb)) .gt. fuzz*10 ,
C     .    dabs(r(3,3)-cb) .gt. fuzz*10
      if (dabs(r(1,1)-(ca*cb*cg-sa*sg)) .gt. fuzz*10 .or.
     .dabs(r(1,2)-(sa*cb*cg+ca*sg)) .gt. fuzz*10 .or.
     .dabs(r(1,3)-(-sb*cg)) .gt. fuzz*10 .or.
     .dabs(r(2,1)-(-ca*cb*sg-sa*cg)) .gt. fuzz*10 .or.
     .dabs(r(2,2)-(-sa*cb*sg+ca*cg)) .gt. fuzz*10 .or.
     .dabs(r(2,3)-(sb*sg)) .gt. fuzz*10 .or.
     .dabs(r(3,1)-(ca*sb)) .gt. fuzz*10 .or.
     .dabs(r(3,2)-(sa*sb)) .gt. fuzz*10 .or.
     .dabs(r(3,3)-cb) .gt. fuzz*10) then
        call awrit3(' rm2eua found alpha=%1;9d  beta=%1;9d'//
     .  '  gamma=%1;9d from rot:',' ',80,i1mach(2),alpha,beta,gamma)
        print 335, ((r(i,j),j=1,3),i=1,3)
  335   format((3f20.15))
        call rx('rm2eua: matrix not consistent with Euler angles: numerical error-->fuzz too large?')
      endif
      end
      subroutine rotma(phi,theta,gamma,rotm)
C- Generate rotation matrix for rotation about specified polar angle
C  Rotates coordinates by gamma about (phi,theta)
C     implicit none
      double precision rotm(3,3),r2(3,3),r1(3,3),phi,theta,gamma,xx
C     integer i,j

      call dpzero(r1,9)
      r1(1,1) = dcos(phi)
      r1(2,2) = r1(1,1)
      r1(2,1) = -dsin(phi)
      r1(1,2) = -r1(2,1)
      r1(3,3) = 1
      call dpzero(r2,9)
      r2(2,2) = 1
      r2(1,1) = dcos(theta)
      r2(3,3) = r2(1,1)
      r2(3,1) = dsin(theta)
      r2(1,3) = -r2(3,1)
c ... rotm <- Rtheta*Rphi
      call dmpy(r2,3,1,r1,3,1,rotm,3,1,3,3,3)
      call dpzero(r2,9)
      r2(1,1) = dcos(gamma)
      r2(2,2) = r2(1,1)
      r2(2,1) = -dsin(gamma)
      r2(1,2) = -r2(2,1)
      r2(3,3) = 1
c ... r1 <- Rgamma*(Rtheta*Rphi)
      call dmpy(r2,3,1,rotm,3,1,r1,3,1,3,3,3)
C ... r2 <- (Rtheta*Rphi)^-1, rotm <- (Rtheta*Rphi)^-1*Rg*(Rtheta*Rphi)
      call dinv33(rotm,0,r2,xx)
      call dmpy(r2,3,1,r1,3,1,rotm,3,1,3,3,3)
C      print 336, ((rotm(i,j),j=1,3),i=1,3)
C  336 format(' rot matrix (Rtheta*Rphi)^-1*Rg*(Rtheta*Rphi):'/(3f12.6))
      end
C      subroutine euler(x,y,z,om,al,be,ga)
C note this convention rotates backwards from usual one.
C      implicit none
C      double precision x,y,z,om,al,be,ga
C      double precision pi,fi,st,ct,ay,tg,om2,r,atg
C      data pi /3.141592653589793d0/
C      om2 = om
C      r = dsqrt(x*x + y*y + z*z)
C      x = x/r
C      y = y/r
C      z = z/r
C      ay = dabs(y)
C      ct = z
C      st = dsqrt(x*x + y*y)
C      if (st .eq. 0d0) then
C        fi = 0d0
C      else
C        fi = dasin(ay/st)
C      endif
C      if (x .gt. 0d0 .and. y .lt. 0d0) fi = 2d0*pi - fi
C      if (x .lt. 0d0 .and. y .gt. 0d0) fi =      pi - fi
C      if (x .lt. 0d0 .and. y .lt. 0d0) fi =      pi + fi
C      be = 2d0*dasin(st*dsin(om2))
C      if (dabs(om/pi-0.5d0) .lt. 0.001d0) then
C        if (ct .ne. 0d0) then
C          atg = pi*dsign(1d0,(0.5d0-om/pi)*ct)/2d0
C        else
C          atg = fi - pi / 2d0
C        end if
C      else
C        tg = ct*dsin(om2)/dcos(om2)
C        atg = datan(tg)
C      endif
C      al = atg + fi - pi/2d0
C      ga = atg - fi + pi/2d0
C      if (ga .lt. -pi) ga = ga+2*pi
C
C      print *, 'fi,be,gp,=',fi,be,atg
C      tg = dtan(atg)
C
C      print *, 'om=',om,datan(tg/ct)
C      print *, dcos(om)**2*tg**2+sin(be/2)**2-dsin(om)**2
C
C      print 1,al,be,ga
C    1 format(' al=',f10.6,' be=',f10.6,' gam=',f10.6)
C      stop
C      end

