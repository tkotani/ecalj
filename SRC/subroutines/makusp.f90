subroutine makusp(n0,z,nsp,rmax,lmxa,v,a,nr,pnu,pnz,rsml,ehl, ul,sl,gz,ruu,rus,rss)
  use m_hansr,only: hansr
!  use m_vxtrap,only: vxtrap
  !- Augmentation functions of pure val,slo (times r) from spherical V and b.c.
  ! ----------------------------------------------------------------------
  !r  ul: linear combination of phi,phidot val=1 slo=0
  !r  sl: linear combination of phi,phidot val=0 slo=1
  !r      ul and sl are r * u and r * s, respectively.
  !r  gz: gz is made for any l for pnz is nonzero where:
  !r      gz = r* ( phi_z - phi_z(rmax) u - (phi_z)'(rmax) s )
  !r      Here phi_z be the w.f. corresponding to loc. orbital spec'd by pnz.
  !r      By construction, gz/r has value=slope = 0 at rmax.
  !i Inputs
  !i   n0    :leading dimension of pnu and pnz
  !i   z     :nuclear charge
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   rmax  :augmentation radius, in a.u.
  !i   lmxa  :augmentation L-cutoff
  !i   v     :spherical potential (atomsr.f)
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
  !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
  !i   pnz   :boundary conditions for (optional) second p.q.n.
  !o Outputs
  !o   ul    :r * linear combination of wave functions; see Remarks
  !o   sl    :r * linear combination of wave functions; see Remarks
  !o   gz    :r * state with 2nd p.q.n; see Remarks
  !o         :If no pnz is nonzero, gz is not touched
  !o   ruu   :diagonal product ul*ul, including small component
  !o         :ruu = ruu(ir,l+1,1,isp)
  !o         :If pnz for a given l channel is nonzero, ruu(ir,l+1,2,isp)
  !o         :is assigned to   gz*ul
  !o   rus   :diagonal product ul*sl, including small component
  !o         :rus = rus(ir,l+1,1,isp)
  !o         :If pnz for a given l channel is nonzero, rus(ir,l+1,2,isp)
  !o         :is assigned to   gz*sl
  !o   rss   :diagonal product sl*sl, including small component
  !o         :rss = rss(ir,l+1,1,isp)
  !o         :If pnz for a given l channel is nonzero, rss(ir,l+1,2,isp)
  !o         :is assigned to  gz*gz
  !l Local variables
  !l   lpzi  :flags how local orbitals is to be treated in current channel
  !l         :0 no local orbital gz
  !l         :1 value and slope of gz constructed to be zero at rmax
  !l         :  by admixture of phi,phidot
  !l         :2 a smooth Hankel tail is attached (extended local orbital)
  !l         :  and it is included explicitly in the basis.
  !l         :3 a smooth Hankel tail is attached (extended local orbital)
  !l         :  and is coupled to the valence states in an extended atom
  !l         :  approximation.
  !b Bugs
  !b   lpzi=3 is not implemented
  !r
  !u Updates
  !u   12 Aug 04 First implementation of extended local orbitals
  !u             Altered argument list.
  !u   21 May 04 (ATP) pass n0 as argument rather than local dimension
  !u   06 Mar 02 Added code to scale gz (see SCALEGZ)
  !u   21 Aug 01 Extended to local orbitals.  Altered argument list.
  !u   16 May 00 Adapted from nfp makusp, makusr and potpsr.
  ! ----------------------------------------------------------------------
  implicit none
  integer :: lmxa,nr,nsp,n0
  double precision :: a,rmax,z,rsml(0:lmxa),ehl(0:lmxa), &
       v(nr,1),pnu(n0,2),pnz(n0,2), &
       ul(nr,lmxa+1,1),sl(nr,lmxa+1,1),gz(nr,lmxa+1,1), &
       ruu(nr,lmxa+1,2,1),rus(nr,lmxa+1,2,1),rss(nr,lmxa+1,2,1)
  integer :: i,l,k,lpz,lpzi,nrbig,idx
!  parameter (nrx=1501)
  double precision :: dphi,dphip,enu,ez,p,phi,phip, &
       g(nr,2),gp(nr,8),gzl(nr,2),phz,dphz
  double precision :: rofi(nr),rwgt(nr)!,vbig(nrx,2)
  double precision :: xi(0:n0),wk(2),fac1
!  logical:: isanrg, l_dummy_isanrg
! --- Make rofi,rwgt, and possibly extended mesh ---
!  call vxtrap(1,z,rmax,lmxa,v,a,nr,nsp,pnz,rs3,vmtz,nrx,lpz,nrbig,rofi,rwgt,vbig)
!  if(nrbig<nr.or.nrbig>nrx) call rx('makusp:nrbig<nr.or.nrbig>nrx')
  call radmsh(rmax,a,nr,rofi)
  do  i = 1, nsp
     do  l = 0, lmxa
        k = l+1
        lpzi = 0
        if (pnz(k,i) >  0)  lpzi = 1
        if (pnz(k,i) >= 10) lpzi = 2
        if (lpzi /= 0) then
           call makrwf(0,z,rmax,l,v(1,i),a,nr,rofi,pnz(1,i),2, gzl,gp,ez,phz,dphz,phip,dphip,p)
           if (lpzi > 1) then !  Scale extended local orbital
              call hansr(rsml(l),0,l,1,[l],[ehl(l)],[rmax**2],1,1,[idx],11,xi)
              fac1 = gzl(nr,1)/rmax/xi(l)
              gzl=gzl/fac1 !call dscal(2*nr,1/fac1,gzl,1)
           endif
        endif
        call makrwf(0,z,rmax,l,v(1,i),a,nr,rofi,pnu(1,i),2,g,gp, &
             enu,phi,dphi,phip,dphip,p)
        !   ... Scale gz so that <|gz-P(g,gp)|^2> = 1
        call makus2(lpzi,nr,rofi,g,gp,gzl,phi,dphi,phip,dphip, &
             phz,dphz,l,enu,ez,z,v(1,i),ul(1,k,i),sl(1,k,i), &
             ruu(1,k,1,i),rus(1,k,1,i),rss(1,k,1,i), &
             ruu(1,k,2,i),rus(1,k,2,i),rss(1,k,2,i))
        if (pnz(k,i) > 0) gz(:,k,i)=gzl(1:nr,1) !call dcopy(nr,gzl,1,gz(1,k,i),1)
     enddo
  enddo
  do  i = 1, nsp
     do  l = 0, lmxa
        if(pnz(l+1,i) == 0) gz(1:nr,l+1,i)=0d0 !zero out missing ones 
     enddo
  enddo
end subroutine makusp
subroutine makus2(lpz,nr,rofi,g,gp,gz,phi,dphi,phip,dphip,phz, &
     dphz,l,e,ez,z,v,  ul,sl,ruu,rus,rss,ruz,rsz,rzz)
  !- Kernel to make u and s from g and gp, for one l
  ! ----------------------------------------------------------------------
  !i Inputs
  !l   lpz   :flags how local orbitals is to be treated in current channel
  !l         :0 no local orbital gz
  !l         :1 value and slope of gz constructed to be zero at rmax
  !l         :  by admixture of phi,phidot
  !l         :2 a smooth Hankel tail is attached (extended local orbital)
  !l         :  but no alteration of gz is made
  !l         :3 same as lpz=2 for this routine.
  !i   nr    :number of radial mesh points
  !i   rofi  :radial mesh points
  !i   g     :normalized wave function times r
  !i   gp    :energy derivative of g
  !i   gz    :(lpz=T) normalized semicore wave function times r
  !i   phi   :wave function at rmax, i.e. g/rmax
  !i   dphi  :slope of wave function at rmax, i.e. d(g/rmax)/dr
  !i   phip  :energy derivative of wave function at rmax
  !i   dphip :energy derivative of slope at rmax
  !i   phz   :value of sc w.f. at rmax: gz/rmax
  !i   dphz  :slope of sc wave function at rmax, i.e. d(gz/rmax)/dr
  !i   l     :l quantum number
  !i   e     :energy eigenvalue, needed for rel. small component
  !i   ez    :energy eigenvalue  of semicore state
  !i   z     :nuclear charge
  !i   v     :spherical potential
  !o Outputs
  !o   ul    :r * linear combination of wave functions; see Remarks
  !o   sl    :r * linear combination of wave functions; see Remarks
  !o   gz    :(lpz=T) returned as (input gz - gz(rmax) u - gz'(rmax) s)
  !o   ruu   :diagonal product u_l*u_l,  including small component
  !o   rus   :diagonal product u_l*s_l,  including small component
  !o   rss   :diagonal product s_l*s_l,  including small component
  !o   ...   The following are made when lpz=T:
  !o   ruz   :diagonal product u_l*gz_l, including small component
  !o   rsz   :diagonal product u_l*gz_l, including small component
  !o   rzz   :diagonal product s_l*gz_l, including small component
  !r Remarks
  !r   This routine makes linear combinations (u,s) of out of phi,phidot
  !r   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
  !r   ul and sl are as r * u and r * s, respectively.
  !u Updates
  !u   13 Jul 04 First implementation of extended local orbitals
  !u   21 Aug 01 extended to computation of semicore states
  ! ----------------------------------------------------------------------
  implicit none
  integer :: l,nr,lpz
  double precision :: dphi,dphip,e,ez,phi,phip,phz,dphz,z
  double precision :: g(nr*2),gp(nr*2),gz(nr*2),rofi(nr),v(nr),ul(nr), &
       sl(nr),ruu(nr),rus(nr),rss(nr),ruz(nr),rsz(nr),rzz(nr)
  integer :: ir,jr
  double precision :: as,au,bs,bu,det,fllp1,gfac,gf11,gf12,gf22,r,sx, &
       tmcr,ux,c,phzl,dphzl
  common /cc/ c
  fllp1 = l*(l+1)
  det = phi*dphip - dphi*phip
  au = dphip/det
  bu = -dphi/det
  as = -phip/det
  bs = phi/det
  if (lpz > 1) then
     phzl = 0
     dphzl = 0
  else
     phzl = phz
     dphzl = dphz
  endif
  if (z /= 0) then
     !       This branch computes products of (g,gp,gz)
     if (lpz /= 0) then
        do  ir = 1, nr
           jr = ir+nr
           r = rofi(ir)
           tmcr = r*c - (r*v(ir) - 2d0*z - r*e)/c
           gf11 = 1d0 + fllp1/tmcr**2
           tmcr = r*c - (r*v(ir) - 2d0*z - r*ez)/c
           gf22 = 1d0 + fllp1/tmcr**2
           gf12 = (gf11 + gf22)/2
           !           (u,s), large component
           ul(ir) = au*g(ir) + bu*gp(ir)
           sl(ir) = as*g(ir) + bs*gp(ir)
           !           (u,s), small component
           ux = au*g(jr) + bu*gp(jr)
           sx = as*g(jr) + bs*gp(jr)
           !           Subtract (phz ul + dphz sl) from gz
           gz(ir) = gz(ir) - phzl*ul(ir) - dphzl*sl(ir)
           gz(jr) = gz(jr) - phzl*ux     - dphzl*sx

           ruu(ir) = gf11*ul(ir)*ul(ir) + ux*ux
           rus(ir) = gf11*ul(ir)*sl(ir) + ux*sx
           rss(ir) = gf11*sl(ir)*sl(ir) + sx*sx

           ruz(ir) = gf12*ul(ir)*gz(ir) + ux*gz(jr)
           rsz(ir) = gf12*sl(ir)*gz(ir) + sx*gz(jr)
           rzz(ir) = gf22*gz(ir)*gz(ir) + gz(jr)*gz(jr)
        enddo

        !       This branch computes products of (g,gp) only
     else
        do  ir = 1, nr
           jr = ir+nr
           r = rofi(ir)
           tmcr = r*c - (r*v(ir) - 2d0*z - r*e)/c
           gfac = 1d0 + fllp1/tmcr**2
           !           (u,s), large component
           ul(ir) = au*g(ir) + bu*gp(ir)
           sl(ir) = as*g(ir) + bs*gp(ir)
           !           (u,s), small component
           ux = au*g(jr) + bu*gp(jr)
           sx = as*g(jr) + bs*gp(jr)

           ruu(ir) = gfac*ul(ir)*ul(ir) + ux*ux
           rus(ir) = gfac*ul(ir)*sl(ir) + ux*sx
           rss(ir) = gfac*sl(ir)*sl(ir) + sx*sx
        enddo
     endif
     ! --- Treat z=0 nonrelativistically ---
  else
     if (lpz /= 0) then
        do  ir = 1, nr
           r = rofi(ir)
           ul(ir) = au*g(ir) + bu*gp(ir)
           sl(ir) = as*g(ir) + bs*gp(ir)
           !           Subtract (phzl ul + dphzl sl) from gz
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
end subroutine makus2

