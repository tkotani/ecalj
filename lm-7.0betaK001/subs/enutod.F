      subroutine enutod(l,rmax,avw,pp,amom,vdel,idmod,pmin,pnu,qnu,eb)
C- Translate band-center energy and shift moments
C ----------------------------------------------------------------
Ci Inputs
Ci   l,rmax,pp
Ci   avw:  average Wigner Seitz radius
Ci   amom: moments corresponding to potential parameters pp
Ci   vdel: difference in origins of energy scales by which the moments
Ci         and enu's are defined:  moments and enu defined with a
Ci         common origin => vdel is 0
Ci   eb:   (idmod=2 only): prescribed shift in enu
Ci   idmod:0 => pnu, qnu are shifted to center of gravity of band
Ci         1 => pnu is frozen; qnu shifted only by vdel.
Ci         2 => pnu, qnu to make change in enu=eb+vdel
Ci         3 => pnu, qnu floated to enu = C
Co Outputs
Co   qnu are shifted to make moments relative to pnu;
Co   pnu may be shifted depending on idmod
Co   eb is shift in energy corresponding to shift in qnu and pnu.
Cr Remarks
Cr   Principal quantum number = integer part of PNU is retained.
Cr   amom and qnu can point to the same address space.
Cr   Note: should be checked for vdel .ne. 0
Cr
Cr   Relation between E and D as follows:
Cr   omeg2 - omeg1 = -S Phi1 Phi2 (D2 - D1), in general,
Cr   so ... for omeg = e - enu, find D(omeg) as follows:
Cr   omeg - omegp = -S Phip Phi(omeg) (D(omeg) - Dp)
Cr   omeg - omegm = -S Phim Phi(omeg) (D(omeg) - Dm)
Cr   Use these to eliminate Phi(omeg) and recall that
Cr   Phim = (2 delta/S sdivw)^1/2 and Phip = Phim sdivw/(2(2l+1) Gamma)
Cr   and  omegm = C - enu and omegp = omegm - delta/Gamma
Cr
Cr   Transformation of the moments:
Cr   The zeroth "moment" m0, i.e. int. phi^2 dE, is related to the
Cr   input moment q0 by
Cr     (1) q0 = m0 + p m2.
Cr   with q0 equal to the the total charge inside the sphere.
Cr   The other moments q1 and q2 are identical to m1 and m2.
Cr
Cr   When shifting by e = change in enu, the wave functions
Cr   phi and phidot are rotated to phit and phitdot as:
Cr     (2) phi    = (phit - e*phitdot)/n,
Cr         phidot = p*e*phit + phitdot)/n
Cr   where
Cr     (3) n^2 = 1 + p*e**2
Cr   This transformation preserves normalization of phi and phidot
Cr   and orthogonality of phi with phidot.  The charge density is
Cr     (4) n(r) := m0*phi**2 + 2*m1*phi*phid + m2*phid**2;
Cr   and the charge density in the shifted representation is
Cr   obtained from the linear transformation (2).  The coefficients
Cr   q0,q1 and q2 corresponding to phit and phitdot are obtained
Cr   from the coefficients q0,q1,q2 related to phi and phidot as:
Cr
Cr     (5) q0t = q0
Cr                    2                                    2
Cr         q1t =  (- e *p*q1 + 2*e*p*q2 - e*q0 + q1)/(1 + e *p)
Cr                    2         2                         2
Cr         q2t =  (- e *p*q2 + e *q0 - 2*e*q1 + q2)/(1 + e *p)
Cr
Cr   When corrections order p are neglected:
Cr
Cr   q0t = q0,  q1t = q1 - e q0,  q2t = e^2 q0 - 2 e q1 + q2
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer idmod,l
      double precision vdel,rmax,pmin
      double precision amom(3),pp(6),PNU,QNU(3),avw,eb
C ... Local parameters
      integer ipqn,iprint
      double precision pi,omeg,omegm,omegp,phim,phip,
     .                 aplus,amins,dnu,sdivw
C ... Local parameters
      double precision q0,q1,q2,d1mach,pfree,dmin,pfuni
      external iprint

      pi = 4.d0*datan(1.d0)
      ipqn = pnu
      if (mod(idmod,10) .ne. 2 .or. dabs(amom(1)) .le. d1mach(3)) eb = 0

C Check if no charge in this channel
      if (dabs(amom(1)) .le. d1mach(3)) goto 50

C Shift enu to c.g. if idmod is zero; otherwise shift moments to enu
      if (mod(idmod,10) .eq. 0) then
        eb = amom(2)/amom(1)
      elseif (mod(idmod,10) .eq. 3) then
        eb = pp(2)-pp(1)
      else
        eb = eb+vdel
      endif

C Skip estimate of pnu if freeze pnu and also already estimate
      if (pnu-ipqn .gt. .001 .and. mod(idmod,10) .eq. 1) goto 50

      pfree = 0.5d0 - datan(dble(l))/pi
      if (pmin .eq. 1) pmin = pfree
      if (pmin .gt. 0) pmin = ipqn + pmin

   40 continue
C     omega(E) = E(P)-enu = C-enu + [(P^-1-gamma)/delta]^-1
C     omegm = omega(D=-l-1) = omega(P=0) = C-enu
C     omegp = omega(D=l) = omega(P^-1=0) = C-enu-delta/gamma
      sdivw = (rmax/avw)**(l+l+1)
      omeg = eb - vdel
      omegm = pp(2) - pp(1)
      omegp = omegm - pp(3)**2/pp(5)
      phim = dsqrt(2/rmax) * pp(3)/dsqrt(sdivw)
      phip = phim / (2*(2*l+1)*pp(5)/sdivw)
      aplus = (omeg-omegp)/phip
      amins = (omeg-omegm)/phim
      dnu = (amins*l + aplus*(l+1)) / (amins-aplus)
      pnu = .5d0 - datan(dnu)/pi + ipqn
C     Case pnu < pmin:  find eb corresponding to pmin
C     To second order, E(P)=enu+omega(D)
      if (pnu .lt. pmin) then
        dmin = -tan(pi*(pmin-ipqn-.5d0))
C       Inverse potential function
        pfuni = sdivw/(2*(2*l+1))*(dmin-l)/(dmin+l+1)
C       P^-1(E) = delta/(E-C)+gamma -> 
C       E-enu = C-enu + [delta/(P^-1-gamma)]
        eb = pp(2)-pp(1) + pp(3)*pp(3)/(pfuni-pp(5))
        pnu = pmin
      endif

C --- Transform moments ---
   50 continue
      if (amom(3) .lt. 0) then
        if (iprint() .ge. 0)
     .  print *, 'enutod: amom(3)<0; resetting to 0 ...'
        amom(3) = 0
      endif
      q0 = amom(1)
      q1 = amom(2) - eb*amom(1)
      q2 = amom(3) - eb*(amom(2) + q1)
C ... Handle special case quadratic approximation to shifts makes q2<0
      if (q2 .lt. 0) then
        eb = 0
        if (iprint() .ge. 10) print *, 'ENUTOD: Set EB=0 to make q2>0'
C Note: this will never happen if eb was zero in the first place
        goto 40
      endif
      qnu(1) = q0
      qnu(2) = q1
      qnu(3) = q2

c     norm = 1 + eb**2*pp(4)
c     n2 =   1 - eb**2*pp(4)
c     qnu(1) = q0
c     qnu(2) = (2*eb*pp(4)*q2 - eb*q0 + q1*n2)/norm
c     qnu(3) = (eb**2*q0 - 2*eb*q1 + q2*n2)/norm
      end

