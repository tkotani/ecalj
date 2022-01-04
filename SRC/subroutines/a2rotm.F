c$$$      subroutine a2rotm(strn,linv,ipr,rotm)
c$$$C- Generate Euler angles from a sequence of rotations, ascii input
c$$$Ci strn: Sequence of rotation matrices.  Syntax: rot1[,rot2...],
c$$$Ci   where each rotj looks like (x,y,z)angle.  x:, y: or z: may be
c$$$Ci   substituted for (x,y,z).  Example:  strn that makes rotation matrix
c$$$Ci   corresponding to the Euler angles alpha=pi/4 beta=pi/3 gamma=pi/2:
c$$$Ci   (0,0,1)pi/4,(0,1,0)pi/3,(0,0,1)pi/2   or just  z:pi/4,y:pi/3,z:pi/2
c$$$Ci linv: T, return transpose of rotm, corresponding to rotation of
c$$$Ci   positions rather than the coordinates (default).
c$$$Ci ipr: verbosity
c$$$Co rotm: rotation matrix.  Follow with a call to rm2eua to find Euler
c$$$Co   angles.
c$$$C     implicit none
c$$$      logical linv
c$$$      character*(*) strn, rchr*7, ss*80
c$$$      double precision rotm(3,3)
c$$$      double precision alpha,beta,gamma,rotj(3,3),r2(3,3),angle,
c$$$     .r(3),pi,ddot,fuzz,ct,theta,phi
c$$$
c$$$      integer ls,itrm,ich,a2vec,ix(3),i,j,irot,ipr,i1mach
c$$$      parameter (fuzz = 1d-12)
c$$$      data rchr /'(XxYyZz'/
c$$$
c$$$
c$$$      pi = 4*datan(1d0)
c$$$      call dpzero(rotm,9)
c$$$      rotm(1,1) = 1
c$$$      rotm(2,2) = 1
c$$$      rotm(3,3) = 1
c$$$      ss = strn
c$$$      ls = len(ss)
c$$$      irot = 0
c$$$      ich = 0
c$$$C --- Entry point to accumulate new rotation ---
c$$$    5 continue
c$$$      irot = irot+1
c$$$      call chrps2(ss,rchr,len(rchr),ich,ich,itrm)
c$$$      if (itrm .eq. 0) goto 999
c$$$      if (itrm .gt. 1) then
c$$$        call dpzero(r,3)
c$$$        if (ss(ich+2:ich+2) .ne. ':') goto 999
c$$$        itrm = itrm/2
c$$$        if (itrm .gt. 3) call rx('bug in a2rotm')
c$$$   11   r(itrm) = 1
c$$$        ich = ich+2
c$$$      else
c$$$        ich = ich+1
c$$$        if (a2vec(ss,ls,ich,4,',)',2,2,3,ix,r) .ne. 3) goto 999
c$$$      endif
c$$$      if (a2vec(ss,ls,ich,4,', ',2,2,1,ix,angle) .ne. 1) goto 999
c$$$      call dpzero(rotj,9)
c$$$C ... Rotation about r
c$$$      call dscal(3,1/dsqrt(ddot(3,r,1,r,1)),r,1)
c$$$      ct = r(3)
c$$$      phi = 0
c$$$      if ( ct .gt.  1 .and.  ct .lt. 1-fuzz) ct = 1
c$$$      if (-ct .gt.  1 .and. -ct .lt. 1-fuzz) ct = -1
c$$$      if (dabs(ct) .lt. 1) phi = datan2(r(2),r(1))
c$$$      theta = dacos(ct)
c$$$      if (theta .lt. 0) theta = 2*pi - theta
c$$$      call rotma(phi,theta,angle,rotj)
c$$$C ... rotm <- rotj*rotm
c$$$      call dpcopy(rotm,r2,1,9,1d0)
c$$$      call dmpy(rotj,3,1,r2,3,1,rotm,3,1,3,3,3)
c$$$      if (ipr .gt. 40) then
c$$$        call rm2eua(rotm,alpha,beta,gamma)
c$$$        call awrit4(' a2rotm:  rotation %i, alpha = %1;6d  '//
c$$$     .  'beta = %1;6d  gamma = %1;6d',' ',80,i1mach(2),
c$$$     .  irot,alpha,beta,gamma)
c$$$        print 335, ((rotm(i,j),j=1,3),i=1,3)
c$$$  335   format((3f15.9))
c$$$
c$$$C ... code to make beta, phi om from Euler angles
c$$$C        st = dsqrt(1-ct**2)
c$$$C        print *, 'beta=',2d0*dasin(st*dsin(angle))
c$$$C        print *, 'theta,phi=',theta,phi
c$$$C        print *, 'phi=',phi,(alpha-gamma-pi)/2
c$$$C        if (dabs((phi-(alpha-gamma-pi)/2)) .gt. 1d-10) then
c$$$C          print *, 'warning: phi off', need fix this
c$$$C          phi = phi - pi
c$$$C        endif
c$$$C        if (dabs((phi-(alpha-gamma-pi)/2)) .gt. 1d-10)  stop 'phi off'
c$$$C        gp = (alpha+gamma)/2
c$$$C        print *, 'gp=',gp
c$$$C        print *, 'angle=', angle,
c$$$C     .    dasin(dsqrt((datan(gp)**2+dsin(beta/2)**2)/(datan(gp)**2+1)))
c$$$C        tg = dcos(theta)*dtan(angle)
c$$$C        tg = dtan(gp)
c$$$C        print *, dcos(angle)**2*tg**2+dsin(beta/2)**2-dsin(angle)**2,
c$$$C     .    dtan(gp),tg
c$$$C        print *, tg**2+dsin(beta/2)**2-dsin(angle)**2*(1+tg**2)
c$$$C        print *, (tg**2+dsin(beta/2)**2)/(1+tg**2)-dsin(angle)**2
c$$$C        print *, dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2))-dsin(angle)
c$$$C        print *, 'angle=', angle,
c$$$C     .    dasin(dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2)))
c$$$CC ...   sin and cos angle, each undetermined up to a sign
c$$$C        sw = dsqrt((tg**2+dsin(beta/2)**2)/(1+tg**2))
c$$$C        cw = dsqrt(1-sw**2)
c$$$C        print *, 'sin angle=', sin(angle),sw
c$$$C        if (tg .lt. 0) cw = -cw
c$$$CC ...   theta, assuming ??? cw positive
c$$$C        th = datan2(dsin(beta/2),cw*tg)
c$$$C        print *, 'theta=', theta, datan2(dsin(beta/2),cw*tg)
c$$$C        print *, 'theta=', theta, datan2(dsin(beta/2),-cw*tg)
c$$$C        st = dsin(th)
c$$$C        ct = dcos(th)
c$$$C        om = datan2(tg,ct)
c$$$C
c$$$Cc ...   check
c$$$C        print *, 'ck tg', tg,tan(om)*ct
c$$$C        print *, 'ck sb', sin(beta/2),sin(th)*sin(om)
c$$$C
c$$$CC ok for ROT=(.1,.2,.3).9*pi
c$$$CC but backwards for ROT=(.1,.2,-.3)-.9*pi ROT=(.1,-.2,-.3)-.9*pi
c$$$C        stop
c$$$C
c$$$C        call euler(r(1),r(2),r(3),angle,alpha,beta,gamma)
c$$$
c$$$      endif
c$$$
c$$$C --- Pick up another rotation matrix, or quit if done ---
c$$$      if (ix(1) .eq. 1) goto 5
c$$$      if (linv) then
c$$$        call dpcopy(rotm,r2,1,9,1d0)
c$$$        do    i = 1, 3
c$$$          do   j = 1, 3
c$$$             rotm(i,j) = r2(j,i)
c$$$          enddo
c$$$          enddo
c$$$        if (ipr .ge. 50) then
c$$$          call rm2eua(rotm,alpha,beta,gamma)
c$$$          call awrit3(' a2rotm:  inverse, alpha = %1;6d  '//
c$$$     .    'beta = %1;6d  gamma = %1;6d',' ',80,i1mach(2),
c$$$     .    alpha,beta,gamma)
c$$$          print 335, ((rotm(i,j),j=1,3),i=1,3)
c$$$        endif
c$$$      endif
c$$$      return
c$$$
c$$$  999 call rx('a2rotm: rotation must be one of x:#,y:#,z:#,(vec)#')
c$$$      end
c$$$      subroutine eua2rm(alpha,beta,gamma,r)
c$$$C- Generate the rotation matrix corresponding to Euler angles
c$$$C ----------------------------------------------------------------------
c$$$C  Definition: R = Rz(gamma)*Ry(beta)*Rz(alpha)
c$$$Cr This definition has the property that a vector v rotated, i.e.
c$$$Cr    R v = v'
c$$$Cr is rotated to v' = zhat 
c$$$Cr when v is defined by the polar coordinates (alpha,beta)
c$$$Cr    vx = cos(alpha) sin(beta)
c$$$Cr    vy = sin(alpha) sin(beta)
c$$$Cr    vx = cos(beta)
c$$$Cr Note also that NB: r^-1 = r+
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$      double precision r(3,3),alpha,beta,gamma
c$$$C     integer i,j
c$$$      double precision ca,cb,cg,sa,sb,sg
c$$$
c$$$      ca = dcos(alpha)
c$$$      sa = dsin(alpha)
c$$$      cb = dcos(beta)
c$$$      sb = dsin(beta)
c$$$      cg = dcos(gamma)
c$$$      sg = dsin(gamma)
c$$$
c$$$C     Euler angles: R3(=R_z(gamma)) R2(=R_y(beta)) R1(=R_z(alpha))
c$$$      r(1,1) =  cg*cb*ca - sg*sa
c$$$      r(2,1) = -sg*cb*ca - cg*sa
c$$$      r(3,1) =  sb*ca
c$$$      r(1,2) =  cg*cb*sa + sg*ca
c$$$      r(2,2) = -sg*cb*sa + cg*ca
c$$$      r(3,2) =  sb*sa
c$$$      r(1,3) = -cg*sb
c$$$      r(2,3) =  sg*sb
c$$$      r(3,3) =  cb
c$$$
c$$$C      print 335, ((r(i,j),j=1,3),i=1,3)
c$$$C  335 format((3f15.9))
c$$$      end
c$$$      subroutine rm2eua(r,alpha,beta,gamma)
c$$$C- Generate Euler angles from a rotation matrix
c$$$C     implicit none
c$$$      double precision r(3,3),alpha,beta,gamma,fuzz,pi
c$$$      integer i1mach,i,j
c$$$      double precision ca,cb,cg,sa,sb,sg,amg,apg,xx,yy,d1mach
c$$$      parameter (fuzz = 1d-8)
c$$$
c$$$      pi = 4*datan(1d0)
c$$$      cb = r(3,3)
c$$$      if (dabs(cb) .gt. 1 .and.  dabs(cb) .lt. 1+fuzz) cb = sign(1d0,cb)
c$$$C --- Handle case beta is near zero or pi ---
c$$$      if (dabs(dabs(cb)-1) .lt. dsqrt(fuzz)) then
c$$$        xx = r(1,3)**2 + r(2,3)**2
c$$$        yy = r(3,1)**2 + r(3,2)**2
c$$$        beta = datan2(dsqrt((xx+yy)/2),r(3,3))
c$$$      else
c$$$        beta = dacos(cb)
c$$$      endif
c$$$      if (dsqrt(r(3,2)**2 + r(3,1)**2) .lt. fuzz .or.
c$$$     .dsqrt(r(2,3)**2 + r(1,3)**2) .lt. fuzz) then
c$$$        if (cb .gt. 0) then
c$$$          apg = datan2(r(1,2),r(2,2))
c$$$          if (r(3,1)**2+r(3,2)**2 .ne. 0d0) then
c$$$            alpha = datan2(r(3,2),r(3,1))
c$$$            gamma = apg - alpha
c$$$          elseif (r(1,3)**2+r(2,3)**2 .ne. 0d0) then
c$$$            gamma = datan2(r(2,3),-r(1,3))
c$$$            alpha = apg - gamma
c$$$          else
c$$$            gamma = 0d0
c$$$            alpha = apg
c$$$          endif
c$$$        else
c$$$          amg = datan2(-r(1,2),r(2,2))
c$$$          if (r(3,1)**2+r(3,2)**2 .ne. 0d0) then
c$$$            alpha = datan2(r(3,2),r(3,1))
c$$$            if (dsqrt(r(3,1)**2+r(3,2)**2) .lt. d1mach(3)*2) alpha = 0
c$$$            gamma = alpha - amg
c$$$          elseif (r(1,3)**2+r(2,3)**2 .ne. 0d0) then
c$$$            gamma = datan2(r(2,3),-r(1,3))
c$$$            alpha = amg + gamma
c$$$          else
c$$$            gamma = 0
c$$$            alpha = amg
c$$$          endif
c$$$        endif
c$$$      else
c$$$        alpha = datan2(r(3,2),r(3,1))
c$$$        gamma = datan2(r(2,3),-r(1,3))
c$$$      endif
c$$$      if (alpha .gt.  pi) alpha = alpha - 2*pi
c$$$      if (alpha .lt. -pi) alpha = alpha + 2*pi
c$$$      if (gamma .gt.  pi) gamma = gamma - 2*pi
c$$$      if (gamma .lt. -pi) gamma = gamma + 2*pi
c$$$C      call awrit3(' rm2eua:  alpha = %1;6d  beta = %1;6d'//
c$$$C     .  '  gamma = %1;6d',' ',80,i1mach(2),alpha,beta,gamma)
c$$$C ... Check validity of entire rotation matrix
c$$$      ca = dcos(alpha)
c$$$      sa = dsin(alpha)
c$$$      cb = dcos(beta)
c$$$      sb = dsin(beta)
c$$$      cg = dcos(gamma)
c$$$      sg = dsin(gamma)
c$$$C      print *,   dabs(r(1,1)-(ca*cb*cg-sa*sg)) .gt. fuzz*10 ,
c$$$C     .    dabs(r(1,2)-(sa*cb*cg+ca*sg)) .gt. fuzz*10 ,
c$$$C     .    dabs(r(1,3)-(-sb*cg)) .gt. fuzz*10 ,
c$$$C     .    dabs(r(2,1)-(-ca*cb*sg-sa*cg)) .gt. fuzz*10 ,
c$$$C     .    dabs(r(2,2)-(-sa*cb*sg+ca*cg)) .gt. fuzz*10 ,
c$$$C     .    dabs(r(2,3)-(sb*sg)) .gt. fuzz*10 ,
c$$$C     .    dabs(r(3,1)-(ca*sb)) .gt. fuzz*10 ,
c$$$C     .    dabs(r(3,2)-(sa*sb)) .gt. fuzz*10 ,
c$$$C     .    dabs(r(3,3)-cb) .gt. fuzz*10
c$$$      if (dabs(r(1,1)-(ca*cb*cg-sa*sg)) .gt. fuzz*10 .or.
c$$$     .dabs(r(1,2)-(sa*cb*cg+ca*sg)) .gt. fuzz*10 .or.
c$$$     .dabs(r(1,3)-(-sb*cg)) .gt. fuzz*10 .or.
c$$$     .dabs(r(2,1)-(-ca*cb*sg-sa*cg)) .gt. fuzz*10 .or.
c$$$     .dabs(r(2,2)-(-sa*cb*sg+ca*cg)) .gt. fuzz*10 .or.
c$$$     .dabs(r(2,3)-(sb*sg)) .gt. fuzz*10 .or.
c$$$     .dabs(r(3,1)-(ca*sb)) .gt. fuzz*10 .or.
c$$$     .dabs(r(3,2)-(sa*sb)) .gt. fuzz*10 .or.
c$$$     .dabs(r(3,3)-cb) .gt. fuzz*10) then
c$$$        call awrit3(' rm2eua found alpha=%1;9d  beta=%1;9d'//
c$$$     .  '  gamma=%1;9d from rot:',' ',80,i1mach(2),alpha,beta,gamma)
c$$$        print 335, ((r(i,j),j=1,3),i=1,3)
c$$$  335   format((3f20.15))
c$$$        call rx('rm2eua: matrix not consistent with Euler angles: numerical error-->fuzz too large?')
c$$$      endif
c$$$      end
c$$$      subroutine rotma(phi,theta,gamma,rotm)
c$$$C- Generate rotation matrix for rotation about specified polar angle
c$$$C  Rotates coordinates by gamma about (phi,theta)
c$$$C     implicit none
c$$$      double precision rotm(3,3),r2(3,3),r1(3,3),phi,theta,gamma,xx
c$$$C     integer i,j
c$$$
c$$$      call dpzero(r1,9)
c$$$      r1(1,1) = dcos(phi)
c$$$      r1(2,2) = r1(1,1)
c$$$      r1(2,1) = -dsin(phi)
c$$$      r1(1,2) = -r1(2,1)
c$$$      r1(3,3) = 1
c$$$      call dpzero(r2,9)
c$$$      r2(2,2) = 1
c$$$      r2(1,1) = dcos(theta)
c$$$      r2(3,3) = r2(1,1)
c$$$      r2(3,1) = dsin(theta)
c$$$      r2(1,3) = -r2(3,1)
c$$$c ... rotm <- Rtheta*Rphi
c$$$      call dmpy(r2,3,1,r1,3,1,rotm,3,1,3,3,3)
c$$$      call dpzero(r2,9)
c$$$      r2(1,1) = dcos(gamma)
c$$$      r2(2,2) = r2(1,1)
c$$$      r2(2,1) = -dsin(gamma)
c$$$      r2(1,2) = -r2(2,1)
c$$$      r2(3,3) = 1
c$$$c ... r1 <- Rgamma*(Rtheta*Rphi)
c$$$      call dmpy(r2,3,1,rotm,3,1,r1,3,1,3,3,3)
c$$$C ... r2 <- (Rtheta*Rphi)^-1, rotm <- (Rtheta*Rphi)^-1*Rg*(Rtheta*Rphi)
c$$$      call dinv33(rotm,0,r2,xx)
c$$$      call dmpy(r2,3,1,r1,3,1,rotm,3,1,3,3,3)
c$$$C      print 336, ((rotm(i,j),j=1,3),i=1,3)
c$$$C  336 format(' rot matrix (Rtheta*Rphi)^-1*Rg*(Rtheta*Rphi):'/(3f12.6))
c$$$      end
c$$$C      subroutine euler(x,y,z,om,al,be,ga)
c$$$C note this convention rotates backwards from usual one.
c$$$C      implicit none
c$$$C      double precision x,y,z,om,al,be,ga
c$$$C      double precision pi,fi,st,ct,ay,tg,om2,r,atg
c$$$C      data pi /3.141592653589793d0/
c$$$C      om2 = om
c$$$C      r = dsqrt(x*x + y*y + z*z)
c$$$C      x = x/r
c$$$C      y = y/r
c$$$C      z = z/r
c$$$C      ay = dabs(y)
c$$$C      ct = z
c$$$C      st = dsqrt(x*x + y*y)
c$$$C      if (st .eq. 0d0) then
c$$$C        fi = 0d0
c$$$C      else
c$$$C        fi = dasin(ay/st)
c$$$C      endif
c$$$C      if (x .gt. 0d0 .and. y .lt. 0d0) fi = 2d0*pi - fi
c$$$C      if (x .lt. 0d0 .and. y .gt. 0d0) fi =      pi - fi
c$$$C      if (x .lt. 0d0 .and. y .lt. 0d0) fi =      pi + fi
c$$$C      be = 2d0*dasin(st*dsin(om2))
c$$$C      if (dabs(om/pi-0.5d0) .lt. 0.001d0) then
c$$$C        if (ct .ne. 0d0) then
c$$$C          atg = pi*dsign(1d0,(0.5d0-om/pi)*ct)/2d0
c$$$C        else
c$$$C          atg = fi - pi / 2d0
c$$$C        end if
c$$$C      else
c$$$C        tg = ct*dsin(om2)/dcos(om2)
c$$$C        atg = datan(tg)
c$$$C      endif
c$$$C      al = atg + fi - pi/2d0
c$$$C      ga = atg - fi + pi/2d0
c$$$C      if (ga .lt. -pi) ga = ga+2*pi
c$$$C
c$$$C      print *, 'fi,be,gp,=',fi,be,atg
c$$$C      tg = dtan(atg)
c$$$C
c$$$C      print *, 'om=',om,datan(tg/ct)
c$$$C      print *, dcos(om)**2*tg**2+sin(be/2)**2-dsin(om)**2
c$$$C
c$$$C      print 1,al,be,ga
c$$$C    1 format(' al=',f10.6,' be=',f10.6,' gam=',f10.6)
c$$$C      stop
c$$$C      end
c$$$
