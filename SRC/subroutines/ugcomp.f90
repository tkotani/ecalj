subroutine ugcomp(qmom,gpot0,hpot0,ugg,f)
  use m_struc_def           !Cgetarg
  use m_lgunit,only:stml
  use m_smhankel,only:hhugbl,hgugbl,ggugbl
  use m_lattic,only: rv_a_opos
  use m_lmfinit,only: ispec,sspec=>v_sspec,nbas
  !- Part of the smooth estatic energy from compensating G's alone.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   sspec :struct containing species-specific information
  !i   slat  :struct containing information about the lattice
  !i   qmom  :multipole moments of on-site densities (rhomom.f)
  ! o Inputs/Outputs
  ! o  Let n0  = smooth potential without compensating gaussians
  ! o      n0~ = smooth potential with compensating gaussians
  ! o    phi0  = ves[n0]
  ! o    phi0~ = ves[n0~]
  ! o    g_RL  = gaussian in RL channel
  ! o    h_R   = l=0 sm hankel in RL channel, for core density
  ! o  Then:
  ! o  gpot0 :On input, integrals g_RL * phi0
  ! o        :On output, integrals g_RL * phi0~
  ! o  hpot0 :On input, integrals h_R * phi0
  ! o        :On output, integrals h_R * phi0~
  !o Outputs
  !o   ugg   :electrostatic energy integral [n0~-n0]*[phi0~-phi0]
  !i   f     :contribution to forces is added
  !r Remarks
  !u Updates
  !u   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   22 Apr 00 Adapted from nfp ugcomp
  ! ----------------------------------------------------------------------
  implicit none
  real(8):: qmom(*) , gpot0(*) , f(3,nbas) , hpot0(nbas) , ugg
  integer :: ndim,ndim0,i,ib,ilm1,ilm2,is,iv0,jb,js,jv0,nvl,l1,l2, &
       lfoc1,lfoc2,ll,lmax1,lmax2,m,nlm1,nlm2
  parameter (ndim=49, ndim0=2)
  double precision :: ceh1,ceh2,cof1,cof2,cofg1,cofg2,cofh1,cofh2,fpi, &
       pi,qcorg1,qcorg2,qcorh1,qcorh2,qsc1,qsc2,qm1,qm2,rg1,rg2,rh1, &
       rh2,srfpi,y0,z1,z2
  double precision :: df(0:20),ff(3),tau1(3),tau2(3)
  double complex s(ndim,ndim),ds(ndim,ndim,3),s0(ndim0,ndim0), &
       ds0(ndim0,ndim0,3),wk(ndim0,ndim0),dwk(ndim0,ndim0,3)
  integer :: nlmx,npmx,ip,mp,nbmx
  parameter (nlmx=64, npmx=1, nbmx=256)
  double precision :: xf(3,nbas,npmx),xhpot0(nbas,npmx), xgpot0(nlmx*nbas,npmx),xugg(npmx)
  integer:: ibini,ibend
  call tcn('ugcomp')
  call stdfac(20,df)
  pi = 4d0*datan(1d0)
  fpi = 4d0*pi
  srfpi = dsqrt(fpi)
  y0 = 1d0/srfpi
  mp = 1
  if (npmx < mp) call rxi('ugcomp: increase npmx, needed',mp)
  ! --- Loop over sites where charge lump making pot is centered ---
  ugg = 0d0
  iv0 = 0
  ip = 1
  call dpzero(xugg, mp)
  call dpzero(xgpot0, nlmx*nbas*mp)
  call dpzero(xf, 3*nbas*mp)
  call dpzero(xhpot0, nbas*mp)
  ibini=1
  ibend=nbas
  do ib=ibini,ibend
     is=ispec(ib) 
     tau1=rv_a_opos(:,ib) 
     lmax1=sspec(is)%lmxl
     rg1=sspec(is)%rg
     call corprm(sspec,is,qcorg1,qcorh1,qsc1,cofg1,cofh1,ceh1,lfoc1, &
          rh1,z1)
     nlm1 = (lmax1+1)**2
     !   ... Loop over sites where charge lump sees the potential
     if (lmax1 > -1) then
        jv0 = 0
        do  jb = 1, nbas
           js=ispec(jb)
           tau2=rv_a_opos(:,jb) 
           lmax2=sspec(js)%lmxl
           rg2=sspec(js)%rg
           if (lmax2 > -1) then
              call corprm(sspec,js,qcorg2,qcorh2,qsc2,cofg2,cofh2,ceh2, lfoc2,rh2,z2)
              nlm2 = (lmax2+1)**2
              if (nlm1 > ndim) call rxi('ugcomp: ndim < nlm1=',nlm1)
              if (nlm2 > ndim) call rxi('ugcomp: ndim < nlm2=',nlm2)
              call ggugbl(tau1,tau2,rg1,rg2,nlm1,nlm2,ndim,ndim,s,ds) 
              ff(1) = 0d0
              ff(2) = 0d0
              ff(3) = 0d0
              do  ilm1 = 1, nlm1
                 l1 = ll(ilm1)
                 qm1 = qmom(iv0+ilm1)
                 if (ilm1 == 1) qm1 = qm1 + y0*(qcorg1-z1)
                 cof1 = qm1*fpi/df(2*l1+1)
                 do  ilm2 = 1, nlm2
                    l2 = ll(ilm2)
                    qm2 = qmom(jv0+ilm2)
                    if (ilm2 == 1) qm2 = qm2 + y0*(qcorg2-z2)
                    cof2 = qm2*fpi/df(2*l2+1)
                    xugg(ip) = xugg(ip) + cof1*cof2*s(ilm1,ilm2)
                    xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip) &
                         + s(ilm1,ilm2)*cof1*fpi/df(2*l2+1)
                    !         ... Forces
                    ff(1) = ff(1) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,1)
                    ff(2) = ff(2) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,2)
                    ff(3) = ff(3) + 0.5d0*cof1*cof2*ds(ilm1,ilm2,3)
                 enddo
              enddo
              !     --- Additional h*h, h*g, g*h terms for foca ---
              if (lfoc1 > 0 .OR. lfoc2 > 0) then
                 call hhugbl(0,tau1,tau2,[rh1],[rh2],[ceh1],[ceh2],1,1,ndim0,ndim0, &
                      wk,dwk,s0,ds0) !slat,
                 xugg(ip) = xugg(ip) + cofh1*s0(1,1)*cofh2
                 xhpot0(jb,ip) = xhpot0(jb,ip) + cofh1*s0(1,1)
                 ff(1) = ff(1) + 0.5d0*cofh1*cofh2*ds0(1,1,1)
                 ff(2) = ff(2) + 0.5d0*cofh1*cofh2*ds0(1,1,2)
                 ff(3) = ff(3) + 0.5d0*cofh1*cofh2*ds0(1,1,3)

                 call hgugbl(tau1,tau2,rh1,rg2,ceh1,1,nlm2,ndim,ndim, &
                      s,ds) !slat,
                 do  ilm2 = 1, nlm2
                    l2 = ll(ilm2)
                    qm2 = qmom(jv0+ilm2)
                    if (ilm2 == 1) qm2 = qm2 + y0*(qcorg2-z2)
                    cof2 = qm2*fpi/df(2*l2+1)
                    xugg(ip) = xugg(ip) + cofh1*s(1,ilm2)*cof2
                    ff(1) = ff(1) + 0.5d0*cofh1*cof2*ds(1,ilm2,1)
                    ff(2) = ff(2) + 0.5d0*cofh1*cof2*ds(1,ilm2,2)
                    ff(3) = ff(3) + 0.5d0*cofh1*cof2*ds(1,ilm2,3)
                    xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip) &
                         + s(1,ilm2)*cofh1*fpi/df(2*l2+1)
                 enddo
                 call hgugbl(tau2,tau1,rh2,rg1,ceh2,1,nlm1,ndim,ndim, s,ds) 
                 do  ilm1 = 1, nlm1
                    l1 = ll(ilm1)
                    qm1 = qmom(iv0+ilm1)
                    if (ilm1 == 1) qm1 = qm1 + y0*(qcorg1-z1)
                    cof1 = qm1*fpi/df(2*l1+1)
                    xugg(ip) = xugg(ip) + cof1*s(1,ilm1)*cofh2
                    ff(1) = ff(1) - 0.5d0*cof1*cofh2*ds(1,ilm1,1)
                    ff(2) = ff(2) - 0.5d0*cof1*cofh2*ds(1,ilm1,2)
                    ff(3) = ff(3) - 0.5d0*cof1*cofh2*ds(1,ilm1,3)
                    xhpot0(jb,ip) = xhpot0(jb,ip) + cof1*s(1,ilm1)
                 enddo
              endif
              if (jb /= ib) then
                 do  m = 1, 3
                    xf(m,ib,ip) = xf(m,ib,ip) - ff(m)
                    xf(m,jb,ip) = xf(m,jb,ip) + ff(m)
                 enddo
              endif
              jv0 = jv0+nlm2
           endif
        enddo
        iv0 = iv0+nlm1
     endif
  enddo
  nvl = iv0
  do  80  ip = 1, mp
     do  82  ib = 1, nbas
        f(1,ib) = f(1,ib) + xf(1,ib,ip)
        f(2,ib) = f(2,ib) + xf(2,ib,ip)
        f(3,ib) = f(3,ib) + xf(3,ib,ip)
        hpot0(ib) = hpot0(ib) + xhpot0(ib,ip)
82   enddo
     do  84  i = 1, nvl
        gpot0(i) = gpot0(i) + xgpot0(i,ip)
84   enddo
     ugg = ugg + xugg(ip)
80 enddo
  call tcx('ugcomp')
end subroutine ugcomp
