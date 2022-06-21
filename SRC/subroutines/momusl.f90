subroutine momusl(z,rmt,lmxa,pnu,pnz,rsml,ehl,lmxl,nlml,a,nr,nsp, rofi,rwgt,v0,v1,qum,vum)
  !- Moments of ul*ul, ul*sl, sl*sl and their integrals with true pot.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   z     :nuclear charge
  !i   rmt   :augmentation radius, in a.u.
  !i   lmxa  :augmentation l-cutoff
  !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
  !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
  !i   lmxl  :l-cutoff for density, potential on the radial mesh
  !i   nlml  :(lmxl+1)*(lmxl+1)
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   nsp   :number of spin channels
  !i   rofi  :radial mesh points
  !i   rwgt  :radial mesh weights: int f(r) dr = sum_ir f(ir) * rwgt(ir)
  !i   v0    :spherical potential defining wave functions, without 2Z/r
  !i   v1    :true potential.  Spherical part need not be v0.
  !o Outputs
  !o   qum   :moments (u,s,gz)_l1 * (u,s,gz)_l2 * r**l
  !o         :qum(l1,l2,l,1) = (u_l1 u_l2) * r**l
  !o         :qum(l1,l2,l,2) = (u_l1 s_l2) * r**l
  !o         :qum(l1,l2,l,3) = (s_l1 s_l2) * r**l
  !o         :qum(l1,l2,l,4) = (u_l1 g_l2) * r**l
  !o         :qum(l1,l2,l,5) = (s_l1 g_l2) * r**l
  !o         :qum(l1,l2,l,6) = (g_l1 g_l2) * r**l
  !o   vum   :integrals ((u,s,gz) * v1 * (u,s,gz)),
  !o         :decomposed by M where v1 = sum_M v1_M Y_M
  !o         :vum(l1,l2,M,1) = (u_l1 v1_M u_l2)
  !o         :vum(l1,l2,M,2) = (u_l1 v1_M s_l2)
  !o         :vum(l1,l2,M,3) = (s_l1 v1_M s_l2)
  !o         :vum(l1,l2,M,4) = (u_l1 v1_M g_l2)
  !o         :vum(l1,l2,M,5) = (s_l1 v1_M g_l2)
  !o         :vum(l1,l2,M,6) = (g_l1 v1_M g_l2)
  !o         :Note that
  !o         :vum(l2,l1,M,2) = (s_l1 v1_M u_l2)
  !o         :vum(l2,l1,M,4) = (g_l1 v1_M u_l2)
  !o         :vum(l2,l1,M,5) = (g_l1 v1_M s_l2)
  !l Local variables
  !l   ul    :radial wave functions; see makusp for definition
  !l   sl    :radial wave functions; see makusp for definition
  !l   ruu   :diagonal w.f. products; see makusp for definition
  !l   rus   :diagonal w.f. products; see makusp for definition
  !l   rss   :diagonal w.f. products; see makusp for definition
  !r Remarks
  !r   Diagonal integrals (l1=l2, l=0 potential) are treated specially:
  !r   small component of w.f. is explicitly taken into account.
  !r   This routine does not rely on any special properties of (u,s,gz)
  !u Updates
  !u   31 Jul 04 normalization of local orbitals now depends on type.
  !u             Altered argument list.
  !u   21 Aug 01 Added local orbitals.  Altered argument list.
  !u   13 Jun 00 spin polarized
  !u   16 May 00 Adapted from nfp momusl.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: lmxa,lmxl,nlml,nr,nsp,n0
  parameter (n0=10)
  double precision :: rmt,a,z,rofi(nr),rwgt(nr),rsml(n0),ehl(n0), &
       v0(nr,nsp),v1(nr,nlml,nsp),pnu(n0,nsp),pnz(n0,nsp), &
       qum(0:lmxa,0:lmxa,0:lmxl,6,nsp), &
       vum(0:lmxa,0:lmxa,nlml,6,nsp)
  logical :: lpz1,lpz2
  integer :: l1,l2,lm,i,mlm,isp
  double precision :: pi,srfpi,xx, &
       quu,qus,qss,quz,qsz,qzz,vuu,vus,vss,vv,vuz,vsz,vzz, &
       ul(nr,0:lmxa,nsp),ruu(nr,0:lmxa,2,nsp), &
       sl(nr,0:lmxa,nsp),rus(nr,0:lmxa,2,nsp), &
       gz(nr,0:lmxa,nsp),rss(nr,0:lmxa,2,nsp)
  call tcn('momusl')
  pi = 4d0*datan(1d0)
  srfpi = dsqrt(4d0*pi)
  gz=0d0
  call makusp(n0,z,nsp,rmt,lmxa,v0,a,nr,xx,xx,pnu,pnz,rsml,ehl, ul,sl,gz,ruu,rus,rss)
  ! --- Moments ---
  do  isp = 1, nsp
     do  l1 = 0, lmxa
        do  l2 = 0, lmxa
           lpz1 = pnz(l1+1,1) .ne. 0
           lpz2 = pnz(l2+1,1) .ne. 0
           do  lm = 0, lmxl
              quu = 0d0
              qus = 0d0
              qss = 0d0
              quz = 0d0
              qsz = 0d0
              qzz = 0d0
              if (l1 == l2 .AND. lm == 0) then
                 do  i = 2, nr
                    xx = rwgt(i)
                    quu = quu + xx * ruu(i,l1,1,isp)
                    qus = qus + xx * rus(i,l1,1,isp)
                    qss = qss + xx * rss(i,l1,1,isp)
                 enddo
                 if (lpz1) then
                    do  i = 2, nr
                       xx = rwgt(i)
                       quz = quz + xx * ruu(i,l1,2,isp)
                       qsz = qsz + xx * rus(i,l1,2,isp)
                       qzz = qzz + xx * rss(i,l1,2,isp)
                    enddo
                 endif
              else
                 do  i = 2, nr
                    xx = rwgt(i) * rofi(i)**lm
                    quu = quu + xx * ul(i,l1,isp) * ul(i,l2,isp)
                    qus = qus + xx * ul(i,l1,isp) * sl(i,l2,isp)
                    qss = qss + xx * sl(i,l1,isp) * sl(i,l2,isp)
                 enddo
                 if (lpz2) then
                    do  i = 2, nr
                       xx = rwgt(i) * rofi(i)**lm
                       quz = quz + xx * ul(i,l1,isp) * gz(i,l2,isp)
                       qsz = qsz + xx * sl(i,l1,isp) * gz(i,l2,isp)
                       qzz = qzz + xx * gz(i,l1,isp) * gz(i,l2,isp)
                    enddo
                 endif
              endif
              qum(l1,l2,lm,1,isp) = quu
              qum(l1,l2,lm,2,isp) = qus
              qum(l1,l2,lm,3,isp) = qss
              qum(l1,l2,lm,4,isp) = quz
              qum(l1,l2,lm,5,isp) = qsz
              qum(l1,l2,lm,6,isp) = qzz
           enddo
        enddo
     enddo
  enddo
  ! --- Integrals ul*ul*V_M, ul*sl*V_M, sl*sl*V_M ---
  do  isp = 1, nsp
     do  mlm = 1, nlml
        do  l1 = 0, lmxa
           do  l2 = 0, lmxa
              lpz1 = pnz(l1+1,1) .ne. 0
              lpz2 = pnz(l2+1,1) .ne. 0
              vuu = 0d0
              vus = 0d0
              vss = 0d0
              vuz = 0d0
              vsz = 0d0
              vzz = 0d0
              if (l1 == l2 .AND. mlm == 1) then
                 do  i = 2, nr
                    vv = v1(i,mlm,isp) - srfpi*2d0*z/rofi(i)
                    vuu = vuu + (rwgt(i)*vv) * ruu(i,l1,1,isp)
                    vus = vus + (rwgt(i)*vv) * rus(i,l1,1,isp)
                    vss = vss + (rwgt(i)*vv) * rss(i,l1,1,isp)
                 enddo
                 if (lpz1) then
                    do  i = 2, nr
                       vv = v1(i,mlm,isp) - srfpi*2d0*z/rofi(i)
                       vuz = vuz + (rwgt(i)*vv) * ruu(i,l1,2,isp)
                       vsz = vsz + (rwgt(i)*vv) * rus(i,l1,2,isp)
                       vzz = vzz + (rwgt(i)*vv) * rss(i,l1,2,isp)
                    enddo
                 endif
              else
                 do  i = 2, nr
                    vv = v1(i,mlm,isp)
                    if (mlm == 1) vv = vv - srfpi*2d0*z/rofi(i)
                    vuu = vuu + (rwgt(i)*vv) * ul(i,l1,isp)*ul(i,l2,isp)
                    vus = vus + (rwgt(i)*vv) * ul(i,l1,isp)*sl(i,l2,isp)
                    vss = vss + (rwgt(i)*vv) * sl(i,l1,isp)*sl(i,l2,isp)
                 enddo
                 if (lpz2) then
                    do  i = 2, nr
                       vv = v1(i,mlm,isp)
                       if (mlm == 1) vv = vv - srfpi*2d0*z/rofi(i)
                       vuz = vuz + (rwgt(i)*vv) * ul(i,l1,isp)*gz(i,l2,isp)
                       vsz = vsz + (rwgt(i)*vv) * sl(i,l1,isp)*gz(i,l2,isp)
                       vzz = vzz + (rwgt(i)*vv) * gz(i,l1,isp)*gz(i,l2,isp)
                    enddo
                 endif
              endif
              vum(l1,l2,mlm,1,isp) = vuu
              vum(l1,l2,mlm,2,isp) = vus
              vum(l1,l2,mlm,3,isp) = vss
              vum(l1,l2,mlm,4,isp) = vuz
              vum(l1,l2,mlm,5,isp) = vsz
              vum(l1,l2,mlm,6,isp) = vzz
           enddo
        enddo
     enddo
  enddo
  call tcx('momusl')
end subroutine momusl
