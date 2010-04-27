      subroutine potpar(nl,nsp,lmx,z,rmax,avw,ekap,lso,loptc,lmpol,a,nr,
     .  rofi,v,pnu,idmod,ves,exc,qnu,idu,uh,jh,thrpv,thrpvl,g,gp,pp,
     .  pprel,sop,pmpol,gradrm)

C- Generates potential parameters for given potential
C  ---------------------------------------------------
Ci Inputs:
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   ekap  :muffin-tin zero plus kap2
Ci   lso   :if true, make spin-orbit coupling parms sop
Ci   loptc :if true, make matrix elements <phi grad phi> etc
Ci   lmpol :if true, make multipole moments of phi,phidot
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   pnu   :Determines b.c. for particular l channel; see remarks
Ci   idmod :0 or 1 : linearization energy enu for l channel determined
Ci         :         by input log derivative spec'd by pnu
Ci         :2      : linearization energy enu for l channel determined
Ci         :         from pp(1)-ves; pnu is output; see remarks
Ci         :See also remarks
Ci   ves   :used only when idmod=2 (see remarks)
Ci   exc   :XC energy used to calc. thrpv
Ci   qnu   :energy-weighted moments of the sphere charges
Ci         :used only in the calculation of 3pV
Ci   idu   :idu(l+1)=0 => this l has no nonlocal U matrix
Ci                   1..3 this l has no nonlocal U matrix
Ci                   4    majority and minority C,E shifted by U; see remarks
Ci                   5    majority and minority C,E shifted by U and J; see remarks
Ci         :Note: this routine doesn't implement LDA+U; see remarks
Ci   uh    :U in LDA+U
Ci   jh    :J in LDA+U
Co Outputs:
Co   pp:   :enu (see remarks), c, srdel, p, gamma (aka q), alpha=gamma
Co         :c = center of band = enu + omega-
Co         :sqrt(delta) has proper sign (that of phi-).
Co         :pp(1) enu
Co         :pp(2) calpha
Co         :pp(3) srdel = sqrt(delta) with proper sign (that of phi-).
Co         :pp(4) palpha
Co         :pp(5) gamma, or Q in Varenna notes
Co         :pp(6) alpha, or qbar in Varenna notes
Co   pprel :potential parameters from the Dirac equation
Co         :pprel(1,l,mu,isp,jsp) cgamma(l,mu,isp,jsp) ?
Co         :pprel(2,l,mu,isp,jsp) gamma(l,mu,isp,jsp)
Co         :pprel(3,l,mu,isp,jsp) delta(l,mu,isp,jsp)
Co         :pprel(4,l,mu,isp,jsp) pgamma(l,mu,isp,jsp)
Co   thrpv : 3 PV for sphere
Co   thrpvl: 3 PV for sphere, decomposed by l-channel and spin
Co   g     : wave function; see rseq
Co   gp    : phidot, except for ...
Co   sop   : spin-orbit coupling parameters
Co   pmpol :integral (phi-or-phidot * phi-or-phidot * r**l) :
Co         :matrix elements of w.f. * wf * r**l for multipole moments
Co   gradrm:
Cr Remarks:
Cr   pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Cr   pp's are generated for a specified potential.
Cr   The enu about which the pp's are generated is calculated from
Cr   the input pnu.  Alternatively, by setting idmod=2 for a specified
Cr   l channel the enu is not calculated but taken from pp(1)-ves, in
Cr   which case pnu is not used.
Cr
Cr   Potential parameters generated from (ekap=0):
Cr     omega(-) = -(phi/phidot) (-l-1-dnu)/(-l-1-dnudot)
Cr     omega(+) = -(phi/phidot) (l-dnu)/(l-dnudot)
Cr     phi(+) = phi + omega(+) phidot
Cr     phi(-) = phi + omega(-) phidot
Cr     C  = e + omega(-)
Cr     Delta = phi(-)**2 sdivw rmax/2
Cr     Gamma = (phi(-)/phi(+)) sdivw / (2 (2l+1))
Cr     where sdivw = (rmax/avw)**(l+l+1)
Cr   These conventions differ by those of MSM by the scale factor
Cr   sdivw for delta and gamma.
Cr   For general ekap (see Kanpur notes)
Cr     C = e - W{K,phi} / W{K,phidot}
Cr       = e - phi/phidot * (D{phi}-D{K})/(D{phidot}-D{K})
Cr     sqrt(delta) = - sqrt(avw/2) / W{K,phidot}
Cr                 = - sqrt(avw/2) /wsr / K_l / phidot /(D{phidot}-D{K})
Cr     gamma = W{J,phidot} / W{K,phidot}
Cr           = J_l/K_l * (D{phidot}-D{J})/(D{phidot}-D{K})
Cr   Uses W{a,b} = r*r*(ab' - a'b) = r*a*b*[D{b}-D{a}]
Cr
Cr   LDA+U parameters:  No LDA+U is implemented in this code.
Cr   For now, parameters are just a device to shift pot pars by a constant.
Cr   You can shift (enu,C) of both spins by U (idu=4)
Cr        or shift (enu,C) of spin1 by U, spin2 by J (idu=5)
Cl Local variables
Cl   lpzi  :flags how local orbitals is to be treated in current channel
Cl         :0 no local orbital gz (none implemented for now)
Cr Updates
Cu  21 Dec 05 (wrl) potential shifts to mimic LDA+U
Cu  11 Jul 05 Updated soprm call to render compatible with fp
Cu  18 Jun 04 (A Chantis) relativistic potential parameters
Cu   4 Apr 04 New matrix elements of <phi Bxc phi>
Cu  05 Jan 04 bug fix making SO matrix elements when nl ne lmx+1
Cu  07 Feb 03 When calc. SO matrix elements, also make <phi|phi>
Cu            for inclusion of external magnetic field.
Cu  20 Apr 98 Optical matrix elements adapted from Sergey Rashkeev
Cu            MvS spin-orbit matrix elements adapted from V. Antropov
C  ---------------------------------------------------
C     implicit none
C Passed parameters
      logical lso,lmpol,loptc
      integer nl,lmx,nr,nsp,idmod(0:nl-1),idu(4)
      double precision uh(4),jh(4)
      double precision z,ekap,rmax,avw,a,ves,thrpv,thrpvl(0:nl-1,nsp)
      double precision rofi(1),v(nr,nsp),g(1),gp(nr,4),pp(6,0:nl-1,nsp),
     .                 pnu(0:nl-1,nsp),qnu(3,0:nl-1,nsp),exc(2),
     .                 sop(0:nl-1,nsp,nsp,9),pmpol(nl,nl,2*nl-1,3,nsp),
     .                 gradrm(4,2,nl,2,nsp),pprel(4,0:nl-1,2*nl,2,2)
C Local parameters
      integer i,j,l,nn,nre,n0,stdo
      parameter (n0=10)
      integer lpzi(0:n0)
      double precision eb1,eb2,pi,b,e,val(5),slo(5),sum,sdivw,
     .  phi,dphi,phip,dphip,dlphi,dlphip,phplus,phmins,
     .  el(2),ql(2),e1me2,e1pe2,e1e2,d,dlmins,dlplus,
     .  fi(0:9),gi(0:9),ptmp1(2,2),ptmp2(2,2),ptmp3(2,2),ptmp4(2,2),
     .  ptmp5(2,2),ptmp6(2,2),ptmp7(2,2),ptmp8(2,2)
C     double precision omegam,omegap
      double precision wk,wkdot,wjdot,scl,enu(0:8,2),fac1,fac2
      integer ophi,odphi,owk,odv,owgt

C Heap allocation
      integer w(1)
      common /w/ w
      double precision tol
      parameter (tol=1d-12)
      external rsq1,rseq,phidx

C ... For Dirac equation
      integer mumax,imu,lrel,nglob
      double precision measur,phisq(0:lmx,nsp),
     .  gmt(2,2),fmt(2,2),gmtde(2,2),fmtde(2,2),gsmt(2,2)
      double precision srd(0:lmx,2), srddot(0:lmx,2)
      double precision frd(0:lmx,2), frddot(0:lmx,2)
      double precision srg(0:lmx,2), srgdot(0:lmx,2)
      double precision frg(0:lmx,2), frgdot(0:lmx,2)

      stdo = nglob('stdo')
      lrel = nglob('lrel')
      do  l = 0, min(lmx,3)
C       LDA+U-like shifts in C
        if (idu(l+1) .le. 3) then
        elseif (idu(l+1) .le. 4) then
          call info2(20,0,0,'%10fpotpar l=%i: enu and C shifted by %d',
     .              l,uh(l+1))
C          if (idu(l+1) .le. 3)
C     .      call rx('POTPAR: LDA+U not implemented in ASA')
        elseif (idu(l+1) .eq. 5) then
          call info5(20,0,0,'%10fpotpar l=%i: enu+ and C+ shifted by'//
     .      ' %d;  enu- and C- shifted by %d',l,uh(l+1),jh(l+1),0,0)
        endif
      enddo

      if (lso .or. lmpol .or. loptc) then
        call defdr(ophi,nr*nl*nsp)
        call defdr(odphi,nr*nl*nsp)
        call defdr(owk,nr*4)
        call defdr(odv,nr)
        call defdr(owgt,nr)
        call radwgt(rmax,a,nr,w(owgt))
        call iinit(lpzi,n0+1)
      endif

      thrpv = 0
      eb1 = -20d0
      eb2 =  20d0
      pi = 4d0*datan(1d0)
      b = rmax/(dexp(a*nr - a) - 1d0)
      call bessl2(ekap*rmax**2,0,lmx+1,fi,gi)

      do  80  i = 1, nsp
        do  10  l = 0, lmx
          if (mod(idmod(l),10) .eq. 2) then
            e = pp(1,l,i) - ves
            call rsq1(0,e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
          else
            e = -0.5d0
            nn = int(pnu(l,i))-l-1
            val(1) = rmax
            slo(1) = 1 + dtan(pi*(0.5d0 - pnu(l,i)))
          endif
          call rseq(eb1,eb2,e,tol,z,l,nn,val,slo,v(1,i),
     .              g,sum,a,b,rofi,nr,nre)
          val(1) = val(1)/dsqrt(sum)
          slo(1) = slo(1)/dsqrt(sum)
          call phidx(1,z,l,v(1,i),0d0,0d0,rofi,nr,2,tol,e,val,slo,nn,
     .      g,gp,phi,dphi,phip,dphip,pp(4,l,i),0d0,0d0,0d0,0d0)
          phisq(l,i) = phi*phi
C     ... Keep local copies of phi and phidot for SO coupling
          if (lso .or. lmpol .or. loptc) then
            call dpscop(g,w(ophi),nr,1,1+nr*(l+nl*(i-1)),1d0)
            call dpscop(gp,w(odphi),nr,1,1+nr*(l+nl*(i-1)),1d0)
            enu(l,i) = e
          endif
          dlphi =  rmax*dphi/phi
          dlphip = rmax*dphip/phip
          sdivw = (rmax/avw)**(l+l+1)
C     ... Following valid for ekap=0 only:
C         omegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
C         omegap = -(phi/phip)*(l-dlphi)/(l-dlphip)
C         phplus = phi + omegap*phip
C         phmins = phi + omegam*phip
C         pp(1,l,i) = e
C         pp(2,l,i) = e + omegam
C         pp(3,l,i) = phmins*dsqrt(sdivw*rmax/2)
C         pp(5,l,i) = phmins/(2*(2*l+1)*phplus)*sdivw
C         pp(6,l,i) = pp(5,l,i)
C     ... The following for general ekap:
C         scl should be 1.  Scale srdel to make compatible w/ ekap=0
          wk     = (dlphi-l)*gi(l)  + (l+l+1)*gi(l+1)
          wkdot  = (dlphip-l)*gi(l) + (l+l+1)*gi(l+1)
          wjdot  = (dlphip-l)*fi(l) + ekap*rmax**2*fi(l+1)/(l+l+1)
C         omegam = -(phi/phip)*wk/wkdot
          scl    = phi*phip*(dlphi-dlphip)*rmax
C     ... Force enu to stay really fixed for idmod=2
          if (mod(idmod(l),10) .ne. 2) then
            pp(1,l,i) = e
          else
            pp(1,l,i) = pp(1,l,i) - ves
          endif
          pp(2,l,i) = e - (phi/phip)*wk/wkdot
          pp(3,l,i) = -dsqrt(sdivw/2/rmax)/phip/wkdot*scl
          pp(5,l,i) = sdivw*wjdot/wkdot
          pp(6,l,i) = pp(5,l,i)
C     ... Constant shift for potential parameters, mimicking LDA+U
          if (idu(l+1) .eq. 4) then
             pp(1,l,i) = pp(1,l,i) + uh(l+1)
             pp(2,l,i) = pp(2,l,i) + uh(l+1)
          endif
C     ... Spin-dependent constant shift for enu,C, mimicking LDA+U
          if (idu(l+1) .eq. 5) then
            if (i .eq. 1) then
              pp(1,l,i) = pp(1,l,i) + uh(l+1)
              pp(2,l,i) = pp(2,l,i) + uh(l+1)
            elseif (i .eq. 2) then
              pp(1,l,i) = pp(1,l,i) + jh(l+1)
              pp(2,l,i) = pp(2,l,i) + jh(l+1)
            endif
          endif
C --- Calculate 3PV ala MSM J. Phys. F 12, 141, Eqn 2.24 ---
C NB: formula valid only for ekap=0
          thrpvl(l,i) = 0
          if (qnu(1,l,i).le.0 .or. qnu(3,l,i).le.0 .or. exc(i) .eq. 0)
     .      goto 20
          d     = qnu(2,l,i)**2 - qnu(1,l,i)*qnu(3,l,i)
          e1pe2 = qnu(2,l,i)*qnu(3,l,i)/d
          e1e2  = qnu(3,l,i)**2/d
C The following two apply if have also third moment ...
C         e1pe2 = (qnu(2,l,i)*qnu(3,l,i) - qnu(1,l,i)*qnu(4,l,i))/d
C         e1e2  = (qnu(3,l,i)**2 - qnu(2,l,i)*qnu(4,l,i))/d
          e1me2 = dsqrt(max(e1pe2**2 - 4*e1e2,0d0))
          if (e1me2 .ne. 0) then
            el(1) = (e1pe2 + e1me2)/2
            el(2) = (e1pe2 - e1me2)/2
            ql(1) = ( qnu(2,l,i) - el(2)*qnu(1,l,i))/e1me2
            ql(2) = (-qnu(2,l,i) + el(1)*qnu(1,l,i))/e1me2
          else
            el(1) = e1pe2
            el(2) = e1pe2
            ql(1) = 0
            ql(2) = 0
          endif

C First order approximation to phi, dphi
          phplus = phi + el(1)*phip
          phmins = phi + el(2)*phip
          dlplus = dlphi - el(1)/phi**2/rmax
          dlmins = dlphi - el(2)/phi**2/rmax

#ifndef FIRST_ORDER_PRESSURE
C Get phi, dphi by reintegrating wave function
          call rsq1(0,el(1)+e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
          dlplus = rmax*(slo(1)-val(1)/rmax)/val(1)
          call gintsr(g,g,a,b,nr,z,el(1)+e,l,v(1,i),rofi,sum)
          phplus = val(1)/rmax/dsqrt(sum)
          call rsq1(0,el(2)+e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
          dlmins = rmax*(slo(1)-val(1)/rmax)/val(1)
          call gintsr(g,g,a,b,nr,z,el(2)+e,l,v(1,i),rofi,sum)
          phmins = val(1)/rmax/dsqrt(sum)
#endif

          thrpvl(l,i) =
     .      rmax*ql(1)*phplus**2*
     .      (dlplus*(dlplus+1) - l*(l+1) + rmax**2*(el(1)+e-exc(i))) +
     .      rmax*ql(2)*phmins**2*
     .      (dlmins*(dlmins+1) - l*(l+1) + rmax**2*(el(2)+e-exc(i)))
          thrpv = thrpv + thrpvl(l,i)

   20     continue


   10   continue
   80 continue

C --- Make spin-orbit matrix elements of phi, phidot ---
      if (lso) then
C       call pshpr(51)

C   ... Gradient of average v
        call dpcopy(v,w(owk),1,nr,1d0/nsp)
        if (nsp .eq. 2) call dpadd(w(owk),v(1,2),1,nr,.5d0)
        call radgra(a,b,nr,rofi,w(owk),w(odv))

C   ... Calculate s-o parameters
        call soprm(5,lpzi,w(ophi),w(odphi),w,nr,nsp,nl-1,lmx,v,w(odv),
     .    enu,w,z,rofi,w(owgt),w(owk),sop,w)

C --- Fully relativistic potential parameters ---
      if (lrel .eq. 2) then
        do  l = 0, lmx
          mumax= 2*(l+1)
          do  imu = 1, mumax
C           mu = dble(imu-l) - 1.5d0
            measur = rmax*0.5d0*(phisq(l,1) + phisq(l,2))
            call rdeq(pp(1,l,1),pp(1,l,2),sop,z,v,rofi,nr,2,a,b,l,lmx,
     .        imu,measur,gmt,fmt,gmtde,fmtde,gsmt,pprel)

           call fdpp(pp(1,l,1),pp(1,l,2),sop,ves,gmt,fmt,gmtde,fmtde,z,
     .        rmax,avw,l,lmx,imu,srd(l,1),srd(l,2),frd(l,1),frd(l,2),
     .        srddot(l,1),srddot(l,2),frddot(l,1),frddot(l,2),
     .        srg(l,1),srg(l,2),frg(l,1),frg(l,2),srgdot(l,1),
     .        srgdot(l,2),frgdot(l,1),frgdot(l,2),pprel,gsmt)

C      ... Test: Give scalar Relativistic small parameter
           if(z .ne. 0) then
             do i = 1, 2
               do j = 1, 2
                 ptmp4(i,j) = 0d0
               enddo
               ptmp4(i,i) = pp(4,l,i)
             enddo
             call clebsh(l,imu,ptmp4,ptmp8)
             do i = 1, 2
               do j = 1, 2
                 pprel(4,l,imu,i,j) = ptmp8(i,j)
               enddo
             enddo
           endif

           if (z .eq. 0) then
             do i = 1, 2
               do j = 1, 2
                 ptmp1(i,j) = 0d0
                 ptmp2(i,j) = 0d0
                 ptmp3(i,j) = 0d0
                 ptmp4(i,j) = 0d0
               enddo
               ptmp1(i,i) = pp(2,l,i)+ves
               ptmp2(i,i) = pp(5,l,i)
               ptmp3(i,i) = pp(3,l,i)
               ptmp4(i,i) = pp(4,l,i)
             enddo
             call clebsh(l,imu,ptmp1,ptmp5)
             call clebsh(l,imu,ptmp2,ptmp6)
             call clebsh(l,imu,ptmp3,ptmp7)
             call clebsh(l,imu,ptmp4,ptmp8)
             do i = 1, 2
               do j = 1, 2
                 pprel(1,l,imu,i,j) = ptmp5(i,j)
                 pprel(2,l,imu,i,j) = ptmp6(i,j)
                 pprel(3,l,imu,i,j) = ptmp7(i,j)
                 pprel(4,l,imu,i,j) = ptmp8(i,j)
C                p2 for the empty spheres
C                pprel(4,l,imu,i,j) = 0.d0
               enddo
             enddo
           endif

          enddo
        enddo
      endif


C   ... Matrix elements for constant magnetic field
        call dvset(w(odv),1,nr,1d0)
        call soprm(2,lpzi,w(ophi),w(odphi),w,nr,nsp,nl-1,lmx,v,w(odv),
     .    enu,w,z,rofi,w(owgt),w(owk),sop(0,1,1,4),w)
C       Correct <dot||dot> term from (already known) diagonal part
C       <phi|phi> should be unity; <phi|dot> should be zero
        do  110  l = 0, lmx
          fac1 = pp(4,l,1)/sop(l,1,1,6)
          fac2 = pp(4,l,2)/sop(l,2,2,6)
          sop(l,1,1,6) = sop(l,1,1,6)*sqrt(fac1*fac1)
          sop(l,1,2,6) = sop(l,1,2,6)*sqrt(fac1*fac2)
          sop(l,2,1,6) = sop(l,2,1,6)*sqrt(fac2*fac1)
          sop(l,2,2,6) = sop(l,2,2,6)*sqrt(fac2*fac2)
  110   continue

C   ... Matrix elements of XC field
        call dpcopy(v,w(odv),1,nr,0.5d0)
        call dpadd(w(odv),v(1,nsp),1,nr,-0.5d0)
C       call prrmsh('vxc ',rofi,w(odv),nr,nr,1)
        call soprm(2,lpzi,w(ophi),w(odphi),w,nr,nsp,nl-1,lmx,v,w(odv),
     .    enu,w,z,rofi,w(owgt),w(owk),sop(0,1,1,7),w)

C       debugging: generate M<B>
C        print *, '!!'
C        if (nsp .eq. 2) then
C          wk = 0
C          do  l = 0, lmx
C            ql(1) = qnu(1,l,1) - pp(4,l,1)*qnu(3,l,1)
C            ql(2) = qnu(1,l,2) - pp(4,l,2)*qnu(3,l,2)
C            fac1 =
C     .        ql(1)*sop(l,1,1,7) +
C     .        qnu(2,l,1)*sop(l,1,1,8) +
C     .        qnu(3,l,1)*sop(l,1,1,9)
C            fac2 =
C     .        ql(2)*sop(l,2,2,7) +
C     .        qnu(2,l,2)*sop(l,2,2,8) +
C     .        qnu(3,l,2)*sop(l,2,2,9)
C            wk = wk + fac1-fac2
C            write(stdo,333) l, fac1,fac2,fac1-fac2
C  333       format(' l=',i1,'   q+<B>=',f12.6,'   q-<B>=',f12.6,
C     .             '   M<B>=',f12.6)
C          enddo
C          write(stdo,'('' <Bxc*M>'',f12.6)') wk
C        endif
      endif

C --- Matrix elements of wave function and its gradient ---
      if (loptc) then
        call rgrme(nsp,nl-1,lmx,nr,rofi,w(owgt),w(ophi),w(odphi),gradrm)
      endif

C --- Multipole moments of phi, phidot ---
      if (lmpol) then
C   ... Orthonormalize phi,phidot
        call soprm(4,lpzi,w(ophi),w(odphi),w,nr,nsp,nl-1,lmx,v,w,enu,w,
     .    z,rofi,w(owgt),w(owk),sop,w)
        call mpint(w(ophi),w(odphi),nl-1,lmx,2*lmx,nr,nsp,
     .    rofi,rofi(1+nr),pmpol)
      endif

      if (lso .or. lmpol) call rlse(ophi)
Cgetarg       end
       end subroutine potpar 

