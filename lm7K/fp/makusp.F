      subroutine makusp(n0,z,nsp,rmax,lmxa,v,a,nr,rs3,vmtz,pnu,pnz,
     .rsml,ehl,ul,sl,gz,ruu,rus,rss)
C- Augmentation fcts of pure val,slo (times r) from spherical V and b.c.
C ----------------------------------------------------------------------
Ci Inputs
Ci   n0    :leading dimension of pnu and pnz
Ci   z     :nuclear charge
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rmax  :augmentation radius, in a.u.
Ci   lmxa  :augmentation L-cutoff
Ci   v     :spherical potential (atomsr.f)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rs3   :minimum smoothing radius for extrapolation of MT potential
Ci         :Not used now
Ci   vmtz  :muffin-tin zero: subtracted from V in the fitting procedure.
Ci         :The asymptotic form of V-vmtz is taken to be zero.
Ci         :Not used now
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pnz   :boundary conditions for (optional) second p.q.n.
Co Outputs
Co   ul    :r * linear combination of wave functions; see Remarks
Co   sl    :r * linear combination of wave functions; see Remarks
Co   gz    :r * state with 2nd p.q.n; see Remarks
Co         :If no pnz is nonzero, gz is not touched
Co   ruu   :diagonal product ul*ul, including small component
Co         :ruu = ruu(ir,l+1,1,isp)
Co         :If pnz for a given l channel is nonzero, ruu(ir,l+1,2,isp)
Co         :is assigned to   gz*ul
Co   rus   :diagonal product ul*sl, including small component
Co         :rus = rus(ir,l+1,1,isp)
Co         :If pnz for a given l channel is nonzero, rus(ir,l+1,2,isp)
Co         :is assigned to   gz*sl
Co   rss   :diagonal product sl*sl, including small component
Co         :rss = rss(ir,l+1,1,isp)
Co         :If pnz for a given l channel is nonzero, rss(ir,l+1,2,isp)
Co         :is assigned to  gz*gz
Cl Local variables
Cl   lpzi  :flags how local orbitals is to be treated in current channel
Cl         :0 no local orbital gz
Cl         :1 value and slope of gz constructed to be zero at rmax
Cl         :  by admixture of phi,phidot
Cl         :2 a smooth Hankel tail is attached (extended local orbital)
Cl         :  and it is included explicitly in the basis.
Cl         :3 a smooth Hankel tail is attached (extended local orbital)
Cl         :  and is coupled to the valence states in an extended atom
Cl         :  approximation.
Cb Bugs
Cb   lpzi=3 is not implemented
Cr Remarks
Cr   This routine makes linear combinations (u,s) of phi,phidot
Cr   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
Cr   ul and sl are r * u and r * s, respectively.
Cr
Cr   Construction of local orbital when val,slo=0 at rmax (lpzi=1):
Cr   Let phi_z be the w.f. corresponding to loc. orbital spec'd by pnz.
Cr   gz is made for any l for pnz is nonzero where:
Cr      gz = r* ( phi_z - phi_z(rmax) u - (phi_z)'(rmax) s )
Cr   By construction, gz/r has value, slope = 0 at rmax.
Cr
Cu Updates
Cu   12 Aug 04 First implementation of extended local orbitals
Cu             Altered argument list.
Cu   21 May 04 (ATP) pass n0 as argument rather than local dimension
Cu   06 Mar 02 Added code to scale gz (see SCALEGZ)
Cu   21 Aug 01 Extended to local orbitals.  Altered argument list.
Cu   16 May 00 Adapted from nfp makusp, makusr and potpsr.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmxa,nr,nsp,n0
      double precision a,rmax,z,rs3,vmtz,rsml(0:lmxa),ehl(0:lmxa),
     .v(nr,1),pnu(n0,2),pnz(n0,2),
     .ul(nr,lmxa+1,1),sl(nr,lmxa+1,1),gz(nr,lmxa+1,1),
     .ruu(nr,lmxa+1,2,1),rus(nr,lmxa+1,2,1),rss(nr,lmxa+1,2,1)
C ... Local parameters
      integer i,l,k,nrx,lpz,lpzi,nrbig,idx
      parameter (nrx=1501)
      double precision dphi,dphip,enu,ez,p,phi,phip,
     .g(nrx,2),gp(nrx,8),gzl(nrx,2),phz,dphz
      double precision rofi(nrx),rwgt(nrx),vbig(nrx,2)
      double precision xi(0:n0),wk(2),fac1


C --- Make rofi,rwgt, and possibly extended mesh ---
      call vxtrap(1,z,rmax,lmxa,v,a,nr,nsp,pnz,rs3,vmtz,nrx,
     .lpz,nrbig,rofi,rwgt,vbig)
      call isanrg(nrbig,nr,nrx,'makusp:','nrbig',.true.)

C --- Loop over spins and l ---
      do  i = 1, nsp
        do  l = 0, lmxa
          k = l+1
          lpzi = 0
          if (pnz(k,i) .gt.  0) lpzi = 1
          if (pnz(k,i) .ge. 10) lpzi = 2
          if (pnz(k,i) .ge. 20) lpzi = 3
          if (lpzi .eq. 3) call rxi('makusp: not implemented lpzi=',lpzi)
          if (lpzi .ne. 0) then
            call makrwf(0,z,rmax,l,v(1,i),a,nr,rofi,pnz(1,i),2,
     .      gzl,gp,ez,phz,dphz,phip,dphip,p)
C         Scale extended local orbital
            if (lpzi .gt. 1) then
c              call hansr(rsml(l),0,l,1,l,ehl(l),rmax**2,1,1,idx,wk,11,xi)
              call hansr(rsml(l),0,l,1,l,ehl(l),rmax**2,1,1,idx,11,xi)
              fac1 = gzl(nr,1)/rmax/xi(l)
              call dscal(2*nr,1/fac1,gzl,1)
            endif
          endif
          call makrwf(0,z,rmax,l,v(1,i),a,nr,rofi,pnu(1,i),2,g,gp,
     .    enu,phi,dphi,phip,dphip,p)

C   ... Scale gz so that <|gz-P(g,gp)|^2> = 1
#if SCALEGZ
          if (lpzi) then
            call ortrwf(0,z,l,v(1,i),nr,nr,nr,rofi,rwgt,enu,enu,ez,
     .      g,gp,gzl,D)
            call dscal(nr*2,1/D,gzl,1)
            phz = phz/D
            dphz = dphz/D
          endif
#endif

          call makus2(lpzi,nr,rofi,g,gp,gzl,phi,dphi,phip,dphip,
     .    phz,dphz,l,enu,ez,z,v(1,i),ul(1,k,i),sl(1,k,i),
     .    ruu(1,k,1,i),rus(1,k,1,i),rss(1,k,1,i),
     .    ruu(1,k,2,i),rus(1,k,2,i),rss(1,k,2,i))
          if (pnz(k,i) .gt. 0) call dcopy(nr,gzl,1,gz(1,k,i),1)
        enddo
      enddo

C --- If at least one semicore state, zero out missing ones ---
      if (lpz .ne. 0) then
        do  i = 1, nsp
          do  l = 0, lmxa
            k = l+1
            if (pnz(k,i) .eq. 0) then
              call dpzero(gz(1,k,i),nr)
            endif
          enddo
        enddo
      endif

C     call prrmsh('ul',rofi,ul,nr,nr,1+lmxa)
C      call wrhomt('ul','ul',0,ul,rofi,nr,1+lmxa,nsp)
C      call wrhomt('sl','sl',0,sl,rofi,nr,1+lmxa,nsp)
C      if (lpz .ne. 0) then
C      call wrhomt('gz','gz',0,gz,rofi,nr,1+lmxa,nsp)
C      endif


      end
      subroutine makus2(lpz,nr,rofi,g,gp,gz,phi,dphi,phip,dphip,phz,
     .dphz,l,e,ez,z,v,ul,sl,ruu,rus,rss,ruz,rsz,rzz)
C- Kernel to make u and s from g and gp, for one l
C ----------------------------------------------------------------------
Ci Inputs
Cl   lpz   :flags how local orbitals is to be treated in current channel
Cl         :0 no local orbital gz
Cl         :1 value and slope of gz constructed to be zero at rmax
Cl         :  by admixture of phi,phidot
Cl         :2 a smooth Hankel tail is attached (extended local orbital)
Cl         :  but no alteration of gz is made
Cl         :3 same as lpz=2 for this routine.
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   g     :normalized wave function times r
Ci   gp    :energy derivative of g
Ci   gz    :(lpz=T) normalized semicore wave function times r
Ci   phi   :wave function at rmax, i.e. g/rmax
Ci   dphi  :slope of wave function at rmax, i.e. d(g/rmax)/dr
Ci   phip  :energy derivative of wave function at rmax
Ci   dphip :energy derivative of slope at rmax
Ci   phz   :value of sc w.f. at rmax: gz/rmax
Ci   dphz  :slope of sc wave function at rmax, i.e. d(gz/rmax)/dr
Ci   l     :l quantum number
Ci   e     :energy eigenvalue, needed for rel. small component
Ci   ez    :energy eigenvalue  of semicore state
Ci   z     :nuclear charge
Ci   v     :spherical potential
Co Outputs
Co   ul    :r * linear combination of wave functions; see Remarks
Co   sl    :r * linear combination of wave functions; see Remarks
Co   gz    :(lpz=T) returned as (input gz - gz(rmax) u - gz'(rmax) s)
Co   ruu   :diagonal product u_l*u_l,  including small component
Co   rus   :diagonal product u_l*s_l,  including small component
Co   rss   :diagonal product s_l*s_l,  including small component
Co   ...   The following are made when lpz=T:
Co   ruz   :diagonal product u_l*gz_l, including small component
Co   rsz   :diagonal product u_l*gz_l, including small component
Co   rzz   :diagonal product s_l*gz_l, including small component
Cr Remarks
Cr   This routine makes linear combinations (u,s) of out of phi,phidot
Cr   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
Cr   ul and sl are as r * u and r * s, respectively.
Cu Updates
Cu   13 Jul 04 First implementation of extended local orbitals
Cu   21 Aug 01 extended to computation of semicore states
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,nr,lpz
      double precision dphi,dphip,e,ez,phi,phip,phz,dphz,z
      double precision g(nr*2),gp(nr*2),gz(nr*2),rofi(nr),v(nr),ul(nr),
     .sl(nr),ruu(nr),rus(nr),rss(nr),ruz(nr),rsz(nr),rzz(nr)
C ... Local parameters
      integer ir,jr
      double precision as,au,bs,bu,det,fllp1,gfac,gf11,gf12,gf22,r,sx,
     .tmcr,ux,c,phzl,dphzl
      common /cc/ c

      fllp1 = l*(l+1)
      det = phi*dphip - dphi*phip
      au = dphip/det
      bu = -dphi/det
      as = -phip/det
      bs = phi/det
      if (lpz .gt. 1) then
        phzl = 0
        dphzl = 0
      else
        phzl = phz
        dphzl = dphz
      endif

      if (z .ne. 0) then

C       This branch computes products of (g,gp,gz)
        if (lpz .ne. 0) then
          do  ir = 1, nr
            jr = ir+nr
            r = rofi(ir)
            tmcr = r*c - (r*v(ir) - 2d0*z - r*e)/c
            gf11 = 1d0 + fllp1/tmcr**2
            tmcr = r*c - (r*v(ir) - 2d0*z - r*ez)/c
            gf22 = 1d0 + fllp1/tmcr**2
            gf12 = (gf11 + gf22)/2
C           (u,s), large component
            ul(ir) = au*g(ir) + bu*gp(ir)
            sl(ir) = as*g(ir) + bs*gp(ir)
C           (u,s), small component
            ux = au*g(jr) + bu*gp(jr)
            sx = as*g(jr) + bs*gp(jr)
C           Subtract (phz ul + dphz sl) from gz
            gz(ir) = gz(ir) - phzl*ul(ir) - dphzl*sl(ir)
            gz(jr) = gz(jr) - phzl*ux     - dphzl*sx

            ruu(ir) = gf11*ul(ir)*ul(ir) + ux*ux
            rus(ir) = gf11*ul(ir)*sl(ir) + ux*sx
            rss(ir) = gf11*sl(ir)*sl(ir) + sx*sx

            ruz(ir) = gf12*ul(ir)*gz(ir) + ux*gz(jr)
            rsz(ir) = gf12*sl(ir)*gz(ir) + sx*gz(jr)
            rzz(ir) = gf22*gz(ir)*gz(ir) + gz(jr)*gz(jr)
          enddo

C       This branch computes products of (g,gp) only
        else
          do  ir = 1, nr
            jr = ir+nr
            r = rofi(ir)
            tmcr = r*c - (r*v(ir) - 2d0*z - r*e)/c
            gfac = 1d0 + fllp1/tmcr**2
C           (u,s), large component
            ul(ir) = au*g(ir) + bu*gp(ir)
            sl(ir) = as*g(ir) + bs*gp(ir)
C           (u,s), small component
            ux = au*g(jr) + bu*gp(jr)
            sx = as*g(jr) + bs*gp(jr)

            ruu(ir) = gfac*ul(ir)*ul(ir) + ux*ux
            rus(ir) = gfac*ul(ir)*sl(ir) + ux*sx
            rss(ir) = gfac*sl(ir)*sl(ir) + sx*sx
          enddo

        endif

C --- Treat z=0 nonrelativistically ---
      else
        if (lpz .ne. 0) then
          do  ir = 1, nr
            r = rofi(ir)
            ul(ir) = au*g(ir) + bu*gp(ir)
            sl(ir) = as*g(ir) + bs*gp(ir)
C           Subtract (phzl ul + dphzl sl) from gz
            gz(ir) = gz(ir) - phzl*ul(ir) - dphzl*sl(ir)

            ruu(ir) = ul(ir)*ul(ir)
            rus(ir) = ul(ir)*sl(ir)
            rss(ir) = sl(ir)*sl(ir)

            ruz(ir) = ul(ir)*gz(ir)
            rsz(ir) = sl(ir)*gz(ir)
            rzz(ir) = gz(ir)*gz(ir)
          enddo
        else
          do  ir = 1, nr
            ul(ir) = au*g(ir) + bu*gp(ir)
            sl(ir) = as*g(ir) + bs*gp(ir)

            ruu(ir) = ul(ir)*ul(ir)
            rus(ir) = ul(ir)*sl(ir)
            rss(ir) = sl(ir)*sl(ir)
          enddo
        endif
      endif

C      call prrmsh('ruu',rofi,ruu,nr,nr,1)
C      call prrmsh('rus',rofi,rus,nr,nr,1)
C      call prrmsh('rss',rofi,rss,nr,nr,1)
C      if (lpz) then
C        call prrmsh('ul',rofi,ul,nr,nr,1)
C        call prrmsh('sl',rofi,sl,nr,nr,1)
C        call prrmsh('gz',rofi,gz,nr,nr,2)
C         call prrmsh('ruz',rofi,ruz,nr,nr,1)
C         call prrmsh('rsz',rofi,rsz,nr,nr,1)
C         call prrmsh('rzz',rofi,rzz,nr,nr,1)
C      endif
      end

