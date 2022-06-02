subroutine phvsfp(mode,nsp,lmxa,ppnl,rmax,sab,vab,hab)
  !- Convert potential parameters to hab,vab,sab
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit nonzero -> make sab
  !i         :10s digit nonzero -> make vab
  !i         :100s digit nonzero -> make hab
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   lmxa  :augmentation l-cutoff
  !i   ppnl  :3rd generation potential parameters; see eg potpus.f
  !i   rmax  :augmentation radius, in a.u.
  !o Outputs
  !o   hab   :<u,s | H | u,s> for each pair uu, us, su, ss; see Remarks
  !o   vab   :<u,s | V | u,s> with V = v + vdif
  !o   sab   :<u,s,gz | 1 | u,s,gz>
  !u Updates
  !u   15 Feb 05 A. Chantis modified from pp2hvs to include gz in sab
  !u   28 Aug 04 sab,vab,hab have parameterized first leading dimension
  !u   28 Mar 01 Created by MvS
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,nsp,lmxa,n0,nppn,nab
  parameter (n0=10,nppn=12,nab=9)
  double precision :: rmax,ppnl(nppn,n0,nsp)
  double precision :: sab(nab,n0,nsp),vab(nab,n0,nsp),hab(nab,n0,nsp)
  integer :: mode0,mode1,mode2,i,la,k
  double precision :: phi,phip,dphi,dphip,dlphi,dlphip,det
  double precision :: au,bu,as,bs,s00,s11,szz,s0z,s1z

  mode0 = mod(mode,10)
  mode1 = mod(mode/10,10)
  mode2 = mod(mode/100,10)
  if (mode1 /= 0 .OR. mode2 /= 0) call rx('pp2hvs: not implemented')

  do  i = 1, nsp
     do  la = 0, lmxa
        k = la+1
        dlphi  = ppnl(3,k,i)/rmax
        dlphip = ppnl(4,k,i)/rmax
        phi    = ppnl(5,k,i)
        phip   = ppnl(6,k,i)
        dphi   = phi*dlphi/rmax
        dphip  = dlphip/rmax*phip
        det    = phi*dphip - dphi*phip
        au = dphip/det
        bu = -dphi/det
        as = -phip/det
        bs = phi/det
        s00 = ppnl(2,k,i)
        s11 = ppnl(7,k,i)
        szz = ppnl(8,k,i)
        s0z = ppnl(9,k,i)
        s1z = ppnl(10,k,i)

        !       s10 = 0d0
        !       s01 = 0d0
        !       sab(1,k,i) = au*s00*au + au*s01*bu + bu*s10*au + bu*s11*bu
        !       sab(2,k,i) = au*s00*as + au*s01*bs + bu*s10*as + bu*s11*bs
        !       sab(3,k,i) = as*s00*au + as*s01*bu + bs*s10*au + bs*s11*bu
        !       sab(4,k,i) = as*s00*as + as*s01*bs + bs*s10*as + bs*s11*bs

        if (mode0 == 1) then
           sab(1,k,i) = au*s00*au + bu*s11*bu
           sab(2,k,i) = au*s00*as + bu*s11*bs
           sab(3,k,i) = as*s00*au + bs*s11*bu
           sab(4,k,i) = as*s00*as + bs*s11*bs
           sab(5,k,i) = szz
           sab(6,k,i) = s0z*au + s1z*bu
           sab(7,k,i) = s0z*as + s1z*bs
        endif

     enddo
  enddo

end subroutine phvsfp


