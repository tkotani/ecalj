module m_smvxcm
  public smvxcm
  contains
subroutine smvxcm(sspec,nbas,lfrce,k1,k2,k3,smrho,&
  smpot,smvxc,smvx,smvc,smexc,repsm,repsmx,repsmc,rmusm,rvmusm, &
       rvepsm,focexc,focex,focec,focvxc,f)
  use m_supot,only: rv_a_ogv,iv_a_okv
  use m_struc_def
  use m_lmfinit,only: rv_a_ocy,   lat_alat, lxcf,nsp
  use m_lattic,only: lat_vol
  use m_supot,only: lat_nabc
  use m_supot,only: lat_ng
  use m_lgunit,only:stdo
  use m_xclda,only: evxcp,evxcv
  !- XC potential for smooth mesh density
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   sspec :struct containing species-specific information
  !i     Elts read:
  !i     Stored:
  !i     Passed to: smcorm smvxc4 corprm
  !i   slat  :struct containing information about the lattice
  !i     Elts read: nabc ng ogv okv vol alat ocy
  !i     Stored:
  !i     Passed to: smcorm
  !i   nbas  :size of basis
  !i   lfrce :1, calculate contribution to forces
  !i   k1,k2,k3 dimensions smooth crystal densities, potentials on a mesh
  !i   smrho :smooth valence density on uniform mesh
  !o Outputs
  !o   smvxc :ex-corr  potential of smoothed density + core corrections
  !o   smvx  :exchange potential of smoothed density + core corrections
  !o   smvc  :correlation potential of smoothed density + core corrections
  !o   smexc :ex-corr  energy density of smoothed density + core corrections
  !o   smpot :smooth total potential; smvxc is added to smpot
  !o   repsm :integrated exchange-correlation energy
  !o   repsmx:integrated exchange energy
  !o   repsmc:integrated correlation energy
  !o   rmusm :int (smrho + smcor1) * vxc[rhosm+smcor1]
  !o         :where smcor1 = portion of core treated directly
  !o   rvmusm:int (smrho) * vxc[rhosm+smcor1]
  !o   rvepsm:int (smrho) * exc[rhosm+smcor1]
  !o   focexc:FOCA exchange-correlation energy:
  !o         :int (smcor2) * vxc[rhosm+smcor1]
  !o         :where smcor2 = portion of core treated perturbatively
  !o   focex :exchange part of focexc
  !o   focec :correlation part of focexc
  !o   focvxc:integral of FOCA exchange-correlation potential:
  !o         :int (smcor2) * (smrho) * (dvxc/drho)
  !o   f     :contribution to forces from xc potential
  !l Local variables
  !l  lxcfun :1s digit sets local xc functional
  !l         :  1    Ceperly-Alder
  !l         :  2    Barth-Hedin (ASW fit)
  !l         :  103 PBE-GGA
  !r Remarks
  !r   smoothed core is partition into core1 + core2.  All atoms with
  !r   lfoc1=1 belong to core1; those with lfoc2=1 belong to core2.
  !r  *core1 is included directly into smrho; the nonlinear XC potential
  !r   is computed from vxc[smrho+smcor1].
  !r  *core2 is included perturbatively: its contribution to the vxc
  !r   is computed from the expansion
  !r     vxc[rho + smcor2] = vxc[rho] + smcor2 * (dvxc/drho)
  !r                       = vxc[rho] + dvxc
  !r   The perturbation correction to int (smrho * vxc) is then
  !r     focvxc = int smrho * smcor2 * (dvxc/drho)
  !r   If the perturbation approach is exact,
  !r     (focvxc+rvmusm) -> rvmusm when computed with smcor2=0
  !r   The corresponding XC energy density is
  !r     exc[rho + smcor2] = exc[rho] + smcor2 * (dexc/drho)
  !r                       = exc[rho] + smcor2 * (vxc-exc)/rho
  !r   The perturbation correction to the XC energy is then
  !r     int smcor2 * (vxc-exc) = focexc - int smcor2 exc[smrho+smcor1]
  !u Updates
  !u   21 Apr 09 Handles GGA functionals
  !u   02 Jul 05 skip sites for which cofh=0
  !u   25 Jun 04 return smexc,rvepsm
  !u   14 Jun 02 rhoex and rhoec (T. Miyake)
  !u    8 Feb 02 smvx and smvc (T. Miyake)
  !u   13 Jun 00 spin polarized
  !u    1 May 00 Adapted from nfp sxc_smooth.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nbas,lfrce,k1,k2,k3,ngabc(3),lxcfun
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  real(8):: f(3,nbas) , repsm(2) , repsmx(2) , repsmc(2) , rmusm(2) &
       , rvmusm(2) , rvepsm(2) , focexc(2) , focex(2) , focec(2) , focvxc(2)
  type(s_spec)::sspec(*)
  double complex smrho(k1,k2,k3,*),smpot(k1,k2,k3,*), &
       smvxc(k1,k2,k3,*),smvx(k1,k2,k3,*),smvc(k1,k2,k3,*), &
       smexc(k1,k2,k3)
  integer:: i , k123 , n1 , n2 , n3 , ng , lfoc1 &
       , lfoc2 , iprint , excsan
  complex(8) ,allocatable :: cgh1_zv(:)
  complex(8) ,allocatable :: cgh2_zv(:)
  complex(8) ,allocatable :: dxcv_zv(:)
  double precision :: vol,sum1,sum2,vxcavg(2),x1,x2,alat
  character(180) :: outs
  integer ::iwdummy
  complex(8),allocatable ::smrho_w(:), smcor_w(:)
  integer:: nnn,isp
  real(8):: sss ,srshift,swmin
  real(8),parameter:: minimumrho=1d-14
  logical::enforce_positive_smrho
  logical:: iprx=.true.
!#if (MPI|MPIK)
  include 'mpif.h'
  integer:: procid=0,ier=0
  integer,parameter::master=0
  call mpi_comm_rank(mpi_comm_world,procid,ier)
  iprx=.false.
  if(procid==master) iprx= .TRUE. 
!endif

  !      stdo = lgunit(1)
  ! hangenglob      nsp  = nglob('nsp')
  !      nsp  = globalvariables%nsp
  ! hangenglob      lxcfun = nglob('lxcf')
  !     lxcfun = globalvariables%lxcf
  call tcn('smvxc')
  lxcfun= lxcf
  ngabc = lat_nabc
  ng    = lat_ng
  vol  = lat_vol
  alat = lat_alat
  vol  = lat_vol
  ! ... Sum of foca hankel heads; break into direct and pert. parts
  !      call defcc (osmrho,-k1*k2*k3*nsp)
  allocate(smrho_w(k1*k2*k3*nsp))
  smrho_w=0d0
  allocate(cgh1_zv(ng), cgh2_zv(ng))
  call smcorm(nbas,sspec,ng,rv_a_ogv,cgh1_zv,cgh2_zv, lfoc1,lfoc2)
  ! ... smrho_w = smrho + smoothed core from foca hankel heads
  k123 = 2*k1*k2*k3
  if(lfoc1 == 1) then
     call gvputf ( ng , 1 , iv_a_okv , k1 , k2 , k3 , cgh1_zv , smrho_w )
     call fftz3(smrho_w,n1,n2,n3,k1,k2,k3,1,0,1)
     if (nsp == 2) then
        call dscal(k123,.5d0,smrho_w,1)
        call dpscop(smrho_w,smrho_w,k123,1,1+k123,1d0)
     endif
     !        call dpadd(smrho_w,smrho,1,k123*nsp,1d0)
     call daxpy(k123*nsp,1d0, smrho,1,smrho_w,1)
  else
     call dpcopy(smrho,smrho_w,1,k123*nsp,1d0)
  endif

  if(enforce_positive_smrho()) then
     !!== negative smrho check== This is also similar with what is done in mkpot.
     !! For GGA, we have to supply positive rho. See enforce_positive_smrho section in mkpot.
     nnn=0
     swmin=0d0
     do i=1,k1*k2*k3*nsp
        sss=dreal(smrho_w(i))
        !        if(sss<0d0) then
        if(sss<minimumrho) then !25July 2011 for the case of He with alat =15Ang (even when all positive, we have NaN for vxc).
           nnn=nnn+1
           if(sss<swmin) then
              swmin=sss
           endif
        endif
     enddo
     if(nnn>0) then
        !         write(6,*) 'smvxcm: negative smrho_w number,min(smrho_w)=',nnn,swmin
        if(iprx) write(6,*) 'smvxcm: smrho_w<minimumrho(jun2011) number,min(smrho_w)=',nnn,swmin !25july2011
        srshift = minimumrho + abs(swmin)
        smrho_w = smrho_w + srshift
        if(iprx) write(6,*) 'smvxcm: enforce positive smrho_w. Add srshift=',srshift
     else
        if(iprx) write(6,*) 'smvxcm: all smrho_w is positive'
     endif
  endif

  ! ... Force density strictly positive definite
  !      print *, 'density strictly pos def?'
  !      call smrpos(w(osmrho),k1,k2,k3,n1,n2,n3)
  rvmusm(2) = 0d0
  focexc(2) = 0d0
  focex(2)  = 0d0
  focec(2)  = 0d0
  focvxc(2) = 0d0

  ! --- Direct branch (lfoc2 .eq. 0) ---
  if (lfoc2 == 0) then
     call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho_w,sum1,sum2)
     call smvxc2 ( 0 ,  nsp , lxcfun , vol , n1 , n2 , n3 ,&
          k1 , k2 , k3 , smrho_w , smvxc , smvx , smvc , smexc  &
          , repsm , repsmx , repsmc , rmusm , vxcavg )
     smpot(:,:,:,1:nsp) = smpot(:,:,:,1:nsp)+smvxc(:,:,:,1:nsp)
     do  i = 1, nsp
        call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i), &
             smrho(1,1,1,i),rvmusm(i),x2)
        call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1), &
             smrho(1,1,1,i),rvepsm(i),x2)
     enddo
     focexc(1) = 0d0
     focex(1)  = 0d0
     focec(1)  = 0d0
     focvxc(1) = 0d0
  endif

  ! --- Perturbation branch (lfoc2 .ne. 0) ---
  if (lfoc2 /= 0) then
     call rx('smvxcm: not suppor perturbation branch now')
     !$$$        allocate(dxcv_zv(k1*k2*k3*nsp))
     !$$$        call smvxc2 ( 1 , slat,nsp , lxcfun , vol , n1 , n2 , n3 , k1 , k2
     !$$$     .  , k3 , smrho_w , smvxc , smvx , smvc , smexc , dxcv_zv
     !$$$     .  , repsm , repsmx , repsmc , rmusm , vxcavg )
     !$$$        call dpadd(smpot,smvxc,1,2*k1*k2*k3*nsp,1d0)
     !$$$        do  i = 1, nsp
     !$$$          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),
     !$$$     .    smrho(1,1,1,i),rvmusm(i),x2)
     !$$$          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),
     !$$$     .    smrho(1,1,1,i),rvepsm(i),x2)
     !$$$        enddo
     !$$$        deallocate(smrho_w)
     !$$$C   ... Assemble core tails for linearized treatment on mesh
     !$$$C       w(osmcor) = portion of core treated perturbatively
     !$$$c        osmcor = osmrho
     !$$$        allocate(smcor_w(k1*k2*k3*nsp))
     !$$$Cchp1         call gvputf ( ng , 1 , w ( okv ) , k1 , k2 , k3 , cgh2_zv
     !$$$Cchp1      .  , smcor_w) !w ( osmcor ) )
     !$$$         call gvputf ( ng , 1 , iv_p_okv , k1 , k2 , k3 , cgh2_zv , smcor_w
     !$$$     .   )
     !$$$        call fftz3(smcor_w,n1,n2,n3,k1,k2,k3,1,0,1)
     !$$$        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smcor_w,sum1,sum2)
     !$$$        call dpzero(focexc,2)
     !$$$        call dpzero(focex,2)
     !$$$        call dpzero(focec,2)
     !$$$        do  i = 1, nsp
     !$$$          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvx(1,1,1,i),smcor_w,
     !$$$     .    x1,x2)
     !$$$          focex(i)  = x1/nsp
     !$$$          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvc(1,1,1,i),smcor_w,
     !$$$     .    x1,x2)
     !$$$          focec(i)  = x1/nsp
     !$$$          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),smcor_w,
     !$$$     .    x1,x2)
     !$$$          focexc(i) = x1/nsp
     !$$$C         Add this term to focexc to make focexc=pert corr to rhov*exc
     !$$$C          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),w(osmcor),
     !$$$C     .      x1,x2)
     !$$$        enddo
     !$$$C       Peturbation correction to smvxc
     !$$$        call smvxc3 ( vol , nsp , n1 , n2 , n3 , k1 , k2 , k3 , smrho
     !$$$     .  , smcor_w , dxcv_zv , smvxc , focvxc )
     !$$$        call dpadd(smpot,smvxc,1,2*k1*k2*k3*nsp,1d0)
     !$$$        if (iprint() .ge. 30) then
     !$$$          outs = ' '
     !$$$          call awrit8('%x   foca'//
     !$$$     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N   foca#%2j#'//
     !$$$     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',
     !$$$     .    outs,len(outs),0,
     !$$$     .    focexc(1)+focexc(2),nsp,focexc,focexc(2),
     !$$$     .    focvxc(1)+focvxc(2),nsp,focvxc,focvxc(2))
     !$$$          call awrit1('%a  charge  =%;12,6D',outs,len(outs),-stdo,sum1)
     !$$$        endif
     !$$$        deallocate(smcor_w,dxcv_zv)
  endif
  deallocate(cgh2_zv,cgh1_zv)
  !      call rlse(osmrho)

  ! --- Force from foca sm-head; cgh1 is workspace ---
  if (lfrce /= 0) then
     allocate(cgh1_zv(ng*nsp))
     call dpzero(f,3*nbas)
     if (lfoc1 > 0 .OR. lfoc2 > 0) then
        call fftz3(smvxc,n1,n2,n3,k1,k2,k3,nsp,0,-1)
        call gvgetf(ng, nsp, iv_a_okv , k1 , k2 , k3 , smvxc , cgh1_zv  )
        call smvxc4(nbas, nsp, sspec, alat, vol, rv_a_ocy, ng , rv_a_ogv , cgh1_zv , f )
     endif
     if (allocated(cgh1_zv)) deallocate(cgh1_zv)
  endif
  call tcx('smvxc')
end subroutine smvxcm


subroutine smvxc2(mode,nsp,lxcfun,vol,n1,n2,n3,k1,k2,k3,smrho, &
     smvxc,smvx,smvc,smexc,       rhoeps,rhoex,rhoec,rhomu,vxcavg)
  use m_ftox
  use m_lgunit,only:stdo
!  use m_vxcnlm,only:vxcnlm
  use m_xclda,only: evxcv,evxcp
  !      use m_struc_def, only: s_lat
  !! Not documented well yet.
  !!= Makes smooth part of xc potential smvxc and optionally dsmvxc/drho =
  ! no perturbation branch (dsmvxc)
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i         :0 do not make dsmvxc/drho
  !i         :1 make dsmvxc/drho
  !i         :10s digit
  !i         : 0 calculated LDA for density as is.
  !i         : 2 for any point where rho<0 or rho_isp<0, zero potential
  !i   nsp   :number of spin channels
  !i  lxcfun :switch defining xc functional (evxcv.f)
  !i         :1s digit sets local xc functional
  !i         :  1    Ceperly-Alder
  !i         :  2    Barth-Hedin (ASW fit)
  !i         :  3,4  LD part of PW91 and PBE
  !i         :100s digit sets gradient corrections
  !i         :  0    LSDA
  !i         :  1    Langreth-Mehl
  !i         :  2    PW91
  !i         :  3    PBE
  !i         :  4    PBE with Becke exchange
  !i   vol   :cell volume
  !i   slat  :struct containing information about the lattice
  !i   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
  !i   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
  !i   smrho :smooth density on uniform mesh
  !o Outputs
  !o   smvxc :xc potential of smoothed density (no core contr.)
  !o   smvx  :exchange potential of smoothed density + core corrections
  !o   smvc  :correlation potential of smoothed density + core corrections
  !o   dsmvxc:dvxc/drho (mode=1)
  !o   rhoeps:integrated exchange-correlation energy
  !o   rhoex :integrated exchange energy
  !o   rhoec :integrated correlation energy
  !o   rhomu :integrated exchange-correlation potential
  !o   vxcavg:average xc potential
  !r Remarks
  !r   For perturbation treatment, take numerical derivatives
  !r   df/dr = d/dr (vxc*r**alfa) instead of d/dr vxc because
  !r   f is nearly linear for alpha=2/3.
  !r
  !r   In the spin polarized case, the smooth core density is not
  !r   spin polarized.  Thus to calc. vxc(rho+drho, m+dm) - vxc(rho,m)
  !r   we use dm=0 and thus drho1 = drho2 = drho/2; thus
  !r     dvxc = lim_drho->0  vxc(rho+drho,rho1+drho/2) - vxc(rho,rho1)
  !r
  !u Updates
  !u   20 Nov 09 New 10s digit for mode
  !u   21 Apr 09 Handles GGA functionals
  !u   14 Jun 02 rhoex and rhoec (T. Miyake)
  !u    8 Feb 02 smvx and smvc (T. Miyake)
  !u   12 Jun 00 spin polarized
  !u    1 May 00 Adapted from nfp vxcd_smooth.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: mode,nsp,k1,k2,k3,n1,n2,n3,lxcfun
  double precision :: rhoeps(2),rhoex(2),rhoec(2),rhomu(2), &
       vxcavg(2),vol
  double complex smvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp), &
       smvx(k1,k2,k3,nsp),smvc(k1,k2,k3,nsp), &
       smexc(k1,k2,k3) !,dsmvxc(k1,k2,k3,nsp)
  ! ... Local parameters
  integer :: i,i1,i2,i3,lxcf,nx,iprint,n1x,lxcg,nglob,mode1
  parameter (n1x=512)
  double precision :: alfa,dfdr,dvdr,f,f1,f2,rrho,fac,dmach,rrmin
  double precision :: repnl(2),rmunl(2),vavgnl(2)
  double precision :: vxc2(n1x,2),vxc1(n1x,2), &
       vx1(n1x,2),vc1(n1x,2), &
       exc2(n1x),exc1(n1x), &
       exc1x(n1x),exc1c(n1x)
  double precision :: rho(n1x),rhos(n1x,2)
  character(180) :: outs
  integer:: ixxx
  real(8):: rhomin
  logical :: newmode=.true.
  integer:: nnn,isp
  real(8):: sss ,smmin(2)
  call tcn('smvxc2')
  if (n1 > n1x) call rxi('smvxc2: increase n1x, need',n1)
  !      stdo = lgunit(1)
  lxcf = mod(lxcfun,100)
  lxcg = mod(lxcfun/100,100)
  alfa = 2d0/3d0
  fac = dmach(1)**(1d0/3d0)
  do  i = 1, 2
     rhoeps(i) = 0
     rhoex(i)  = 0
     rhoec(i)  = 0
     rhomu(i)  = 0
     vxcavg(i) = 0
  enddo
  rrmin = 0
  nx = 0
  mode1 = mod(mode/10,10)

  smvxc=0d0
  smvx=0d0
  smvc=0d0
  smexc=0d0
  !     call dpzero(smvxc,2*k1*k2*k3*nsp)
  !     call dpzero(smvx,2*k1*k2*k3*nsp)
  !     call dpzero(smvc,2*k1*k2*k3*nsp)
  !     call dpzero(smexc,2*k1*k2*k3)

  !     Vector of points for each i2,i3.
  do  i3 = 1, n3
     do  i2 = 1, n2

        !     ... rho = total rho for this vec; rhos = spin pol rho
        call dcopy(n1,smrho(1,i2,i3,1),2,rho,1)
        call dcopy(n1,smrho(1,i2,i3,1),2,rhos,1)
        if (nsp == 2) then
           call dcopy(n1,smrho(1,i2,i3,2),2,rhos(1,2),1)
           call daxpy(n1,1d0,rhos(1,2),1,rho,1)
        endif

        !! ... Put df/dr into dsmvxc <--- not support now
        if (mod(mode,10) /= 0) then
           call rx('svxvcm: not support perturbation treatment anymore')
        endif

        !     ... Exchange into smvxc
        if (lxcf == 3 .OR. lxcf == 4) then
           call evxcp(rhos(1,1),rhos(1,2),n1,nsp,lxcf,exc1x,exc1c,exc1, &
                vx1(1,1),vx1(1,2),vc1(1,1),vc1(1,2),vxc1(1,1),vxc1(1,2))
        else
           do  i = 1, nsp
              call evxcv(rho,rhos(1,i),n1,nsp,lxcf,exc1,exc1x,exc1c, &
                   vxc1(1,i),vx1(1,i),vc1(1,i))
           enddo
        endif
        if (mode1 == 2) then
           do  i = 1, nsp
              do  i1 = 1, n1
                 if (rho(i1) <= 0 .OR. rhos(i1,i) <= 0) then
                    vxc1(1,i) = 0
                 endif
              enddo
           enddo
        endif

        do  i = 1, nsp
           !            call evxcv(rho,rhos(1,i),n1,nsp,lxcf,exc1,
           !     .      exc1x,exc1c,vxc1(1,i),vx1(1,i),vc1(1,i))
           call dcopy(n1,vxc1(1,i),1,smvxc(1,i2,i3,i),2)
           call dcopy(n1,vx1(1,i),1,smvx(1,i2,i3,i),2)
           call dcopy(n1,vc1(1,i),1,smvc(1,i2,i3,i),2)
        enddo
        call dcopy(n1,exc1,1,smexc(1,i2,i3),2)

        ! CC... Perturbation dv/dr into dsmvxc <---'svxvcm: not support perturbation treatment anymore'
        !$$$          if (mod(mode,10) .ne. 0) then
        !$$$            do  i = 1, nsp
        !$$$              do  i1 = 1, n1
        !$$$                rrho = rho(i1)
        !$$$                if (rrho .gt. 0) then
        !$$$                  f = vxc1(i1,i) * rrho**alfa
        !$$$                  dvdr = (vxc2(i1,i) - alfa*f/rrho) / rrho**alfa
        !$$$                  dsmvxc(i1,i2,i3,i) = dvdr
        !$$$                else
        !$$$                  dsmvxc(i1,i2,i3,i) = 0
        !$$$                endif
        !$$$              enddo
        !$$$            enddo
        !$$$          endif

        !     ... Add to integrals
        do  i = 1, nsp
           do  i1 = 1, n1
              rrho = rhos(i1,i)
              rrmin = min(rrho,rrmin)
              if (rrho < 0d0) nx = nx+1
              rhomu(i)  = rhomu(i)  + rrho*vxc1(i1,i)
              rhoeps(i) = rhoeps(i) + rrho*exc1(i1)
              rhoex(i)  = rhoex(i)  + rrho*exc1x(i1)
              rhoec(i)  = rhoec(i)  + rrho*exc1c(i1)
              vxcavg(i) = vxcavg(i) + vxc1(i1,i)
           enddo
        enddo
     enddo
  enddo

  f = vol/(n1*n2*n3)
  do  i = 1, nsp
     rhoeps(i) = rhoeps(i)*f
     rhoex(i) = rhoex(i)*f
     rhoec(i) = rhoec(i)*f
     rhomu(i) = rhomu(i)*f
     vxcavg(i) = vxcavg(i)/(n1*n2*n3)
  enddo

  !      do  i = 1, nsp
  !        call zprm3('LDA smvxc (isp=%i)',i,smvxc(1,1,1,i),n1,n2,n3)
  !        call zprm3('dsmvxc (isp=%i)',i,dsmvxc(1,1,1,i),n1,n2,n3)
  !      enddo
  !      print *,' goto vxcnlm lxcg=',lxcg,sum(abs(smrho)),sum(smrho)

  ! ... Gradient corrections to potential, energy
  if(lxcg/=0.and.lxcg/=1) call rx('smvxc2: lxcg/=0.and.lxcg/=1')
  if (lxcg == 1) then
     !        print *,'smvxcm:calling vxcnlm sum smrho=',sum(abs(smrho))
     call vxcnlm(lxcg,nsp,k1,k2,k3,smrho, repnl,rmunl,vavgnl,smvx,smvc,smvxc)
!     newmode=.true.
!     if(newmode) then
        do i=1,nsp
           rhoeps(i) =  repnl(i)
           rhomu(i)  =  rmunl(i)
           vxcavg(i) =  vavgnl(i)
           repnl(i)  = 0d0
           rmunl(i)  = 0d0
           vavgnl(i) = 0d0
        enddo
!     endif
     !        repnl = 0 ; rmunl = 0 ; vavgnl = 0; print *, '!!'
  endif
  !      do  i = 1, nsp
  !        call zprm3('LDA+GGA smvxc (isp=%i)',i,smvxc(1,1,1,i),n1,n2,n3)
  !      enddo

  !     call setpr(30)

  ! ... LDA, GGA Printout
  if (nx > 0) call info5(20,0,0,' smvxcm (warning) mesh density ' &
       //'negative at %i point%?#n>1#s##:  rhomin=%;3g',nx,nx,rrmin,0,0)

  !$$$      if (iprint() .ge. 30 .and. lxcg .ne. 0) then
  !$$$        call awrit8('%x sm GGA'//
  !$$$     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%1f#%2j#sm GGA'//
  !$$$     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
  !$$$     .    0,rhoeps(1)+rhoeps(2),nsp,rhoeps,rhoeps(2),rhomu(1)+rhomu(2),
  !$$$     .    nsp,rhomu,rhomu(2))
  !$$$        call awrit5('%a%?#n==2#%N%3f#  #avg GGA vxc ='//
  !$$$     .    '%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
  !$$$     .    -stdo,nsp,(vxcavg(1)+vxcavg(nsp))/2,nsp,vxcavg,vxcavg(2))
  !$$$        call awrit8('%x sm GGA'//
  !$$$     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%1f#%2j#sm GGA'//
  !$$$     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
  !$$$     .    0,repnl(1)+repnl(2),nsp,repnl,repnl(2),rmunl(1)+rmunl(2),
  !$$$     .    nsp,rmunl,rmunl(2))
  !$$$        call awrit5('%a%?#n==2#%N%3f#  #avg GGA vxc ='//
  !$$$     .    '%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
  !$$$     .    -stdo,nsp,(vavgnl(1)+vavgnl(nsp))/2,nsp,vavgnl,vavgnl(2))
  !$$$      endif
  if (lxcg /= 0) then
     do  i = 1, nsp
        rhoeps(i) = rhoeps(i) + repnl(i)
        !         rhoex(i) = rhoex(i) +
        !         rhoec(i) = rhoec(i) +
        rhomu(i) = rhomu(i) + rmunl(i)
        vxcavg(i) = vxcavg(i) + vavgnl(i)
     enddo
  endif

  ! ... Printout, total potential
  if (iprint() >= 30) then
     !$$$        call awrit8('%x smooth'//
     !$$$     .  ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%7f#%2j#'//
     !$$$     .  '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
     !$$$     .  0,rhoeps(1)+rhoeps(2),nsp,rhoeps,rhoeps(2),rhomu(1)+rhomu(2),
     !$$$     .  nsp,rhomu,rhomu(2))
     !$$$        call awrit5('%a%?#n==2#%N%7f#  #'//
     !$$$     .  'avg vxc =%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
     !$$$  .  -stdo,nsp,(vxcavg(1)+vxcavg(nsp))/2,nsp,vxcavg,vxcavg(2))
     do i=1,nsp
        write(stdo,ftox)'smooth isp rhoeps rhomu vxcavg=' &
             ,i,ftof(rhoeps(i)),ftof(rhomu(i)),ftof(vxcavg(i))
     enddo
  endif

  !      call zprm3('smvxc',0,smvxc,n1,n2,n3)
  !      call zprm3('dsmvxc',0,dsmvxc,n1,n2,n3)
  !      if (nsp .eq. 2) then
  !        call zprm3('smvxc spin 2',0,smvxc(1,1,1,2),n1,n2,n3)
  !        call zprm3('dsmvxc spin 2',0,dsmvxc(1,1,1,2),n1,n2,n3)
  !      endif

  call tcx('smvxc2')
end subroutine smvxc2


subroutine smvxc3(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,smcor,dsmvxc, &
     smvxc,rmuxcc)
  use m_lgunit,only:stdo

  !- Smooth core density times dvxc/drho
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   vol   :cell volume
  !i   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
  !i   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
  !i   smrho :smooth density on n1,n2,n3 mesh
  !i   smcor :smooth core density on n1,n2,n3 mesh
  !i   dsmvxc:dvxc/drho on n1,n2,n3 mesh mesh
  !o Outputs
  !o   smvxc :(dvxc/drho * smcor) = pert. correction to expansion
  !o         : vxc[rho + rhoc] = vxc[rho] + rhoc * dvxc/drho
  !o   rmuxcc:integral smrho * (dvxc/drho * smcor)
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nsp,k1,k2,k3,n1,n2,n3
  double precision :: rmuxcc(nsp),vol
  double complex smvxc(k1,k2,k3,nsp),smcor(k1,k2,k3), &
       dsmvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp)
  ! ... Local parameters
  integer :: i,i1,i2,i3
  double complex cadd,csum(2)

  rmuxcc(2) = 0
  do  i = 1, nsp
     csum(i) = 0d0
     do  i3 = 1, n3
        do  i2 = 1, n2
           do  i1 = 1, n1
              cadd = dsmvxc(i1,i2,i3,i)*smcor(i1,i2,i3)
              smvxc(i1,i2,i3,i) = cadd
              csum(i) = csum(i) + smrho(i1,i2,i3,i)*cadd
           enddo
        enddo
     enddo
     csum(i) = csum(i)*vol/(n1*n2*n3)
     rmuxcc(i) = dble(csum(i))
  enddo

  !     write(stdo,862) csum
  ! 862 format(' csum=',2f14.8)

end subroutine smvxc3

subroutine smvxc4(nbas,nsp,sspec,alat,vol,cy,ng,gv,cvxc,f)
  use m_lgunit,only:stdo
  use m_struc_def  !Cgetarg
  use m_lattic,only: rv_a_opos
  use m_lmfinit,only: ispec
  !- For foca, adds force from shift of smH-head against Vxc.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   nsp   :number of spin channels
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: *
  !i     Stored:    *
  !i     Passed to: corprm
  !i   cy    :Normalization constants for spherical harmonics
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !i   cvxc  :Fourier transform of smooth vxc potential.
  !o Outputs
  !o   f     :force from shift of smH-head against Vxc added to f.
  !r Remarks
  !u Updates
  !u   02 Jul 05  skip sites for which cofh=0
  !u    1 May 00  Adapted from nfp smc_force.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nbas,nsp,ng
  real(8):: gv(ng,3) , alat , vol , cy(1) , f(3,nbas)
  type(s_spec)::sspec(*)

  double complex cvxc(ng,nsp)
  ! ... Local parameters
  integer :: k0,nlmx,kmax,ib,is,lfoc,i,kb,iprint
  double precision :: tau(3),v(3),pi,tpiba,qcorg,qcorh,qsc,cofg, &
       cofh,ceh,rfoc,z,sum1,sum2,sum3,xx
  parameter (k0=3, nlmx = 9)
  double complex gkl(0:k0,nlmx),ccc,cvxci

  !      stdo = lgunit(1)
  pi = 4d0*datan(1d0)
  tpiba = 2d0*pi/alat
  kmax = 0

  ! --- Loop over sites ---
  if (iprint() >= 50) write(stdo,400)
  do  ib = 1, nbas
     is=ispec(ib)
     tau=rv_a_opos(:,ib) 
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     if (lfoc > 0 .AND. cofh /= 0) then
        sum1 = 0d0
        sum2 = 0d0
        sum3 = 0d0
        do  i = 1, ng
           v(1) = gv(i,1)
           v(2) = gv(i,2)
           v(3) = gv(i,3)
           call hklft(v,rfoc,ceh,tau,alat,kmax,1,k0,cy,gkl)
           ccc = cofh*gkl(0,1)/vol
           cvxci = 0.5d0 * (cvxc(i,1) + cvxc(i,nsp))
           xx = -dimag(dconjg(cvxci) * ccc)
           sum1 = sum1 + xx*gv(i,1)
           sum2 = sum2 + xx*gv(i,2)
           sum3 = sum3 + xx*gv(i,3)
        enddo
        sum1 = sum1*vol*tpiba
        sum2 = sum2*vol*tpiba
        sum3 = sum3*vol*tpiba
        f(1,ib) = f(1,ib) + sum1
        f(2,ib) = f(2,ib) + sum2
        f(3,ib) = f(3,ib) + sum3
        do  kb = 1, nbas
           f(1,kb) = f(1,kb) - sum1/nbas
           f(2,kb) = f(2,kb) - sum2/nbas
           f(3,kb) = f(3,kb) - sum3/nbas
        enddo
     endif
  enddo
  if (iprint() >= 50) &
       write(stdo,340) (ib,f(1,ib),f(2,ib),f(3,ib),ib = 1,nbas)
340 format(i4,3f12.6)
400 format(/' xc-force from foca:')

end subroutine smvxc4

!!= Gradient correction to smoothed rho(q) tabulated on a mesh =
!!*Kotani's version newmode with xcpbe.F in abinit Aug2010
subroutine vxcnlm(lxcg,nsp,k1,k2,k3,smrho,repnl,rmunl,vavgnl,vxnl,vcnl,vxcnl)
  use m_supot,only: iv_a_okv,rv_a_ogv
  !      use m_struc_def, only:s_lat
  use m_xcpbe,  only: xcpbe
  use m_lmfinit,only:    lat_alat
  use m_lattic,only: lat_vol
  use m_supot,only: lat_nabc
  use m_supot,only: lat_ng
  !! ----------------------------------------------------------------------
  ! i Inputs
  ! i   lxcg  : dummy now.  (need to set option in xcpbe)
  ! i   slat,smrho(k1,k2,k3,nsp)
  ! o Outputs (for newmode=T).
  ! o   repnl : integral smrho * eps
  ! o   rmunl : integral smrho * vxc
  ! o   vavgnl:average NL XC potential
  ! o   vxcnl : XC potential on uniform mesh.
  !!   vcnl  : dummy (it was correlation part of vxcnl)
  !!   vxnl  : dummy (it was exchange part of vxcnl)
  !! ----------------------------------------------------------------------

  ! cccccccccccccccccccccccccc
  !  old document below. Kink can exist for (grad |grad rho|) (imagine a case with rho=x^2+1)

  ! cccccccccccccSpecifies GGA for old case
  !i         :  0    LSDA
  !i         :  1    Langreth-Mehl
  !i         :  2    PW91
  !i         :  3    PBE
  !i         :  4    PBE with Becke exchange
  !i   nsp   : 2 for spin-polarized case, otherwise 1
  !i   k1..k3: dimensions of smrho,vnl for smooth mesh density
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: nabc ng ogv okv alat vol
  !i     Stored:
  !i     Passed to: vxcgga vxnlcc vxnloc
  !i   smrho :smooth density on uniform mesh
  !l Local variables :
  !l   agr(*,1)  : |grad rhop| or |grad rho| if nsp=1
  !l   agr(*,2)  : |grad rhom| (nsp=2)
  !l   agr(*,k)  : |grad total rho|. k=3 for nsp=2; else k=1
  !l   agr(*,4)  : grad rho+ . grad rho- (only for Langreth-Mehl-Hu)
  !l   ggr(*,1)  : Laplacian of rhop (total rho if nsp=1)
  !l   ggr(*,2)  : Laplacian of rhom (nsp=2)
  !l   gagr(*,k) : (grad rho).(grad |grad rho|)
  !l   gagr(*,1) : (grad rhop).(grad |grad rhop|) (total rho if nsp=1)
  !l   gagr(*,2) : (grad rhom).(grad |grad rhom|) (nsp=2)
  !l   gagr(*,k) : (grad rho).(grad |grad rho|). k=3 for nsp=2; else k=1
  !r Remarks
  !r
  !u Updates
  !u   06 Apr 09 Adapted from vxcnlp.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: lxcg,k1,k2,k3,nsp
  real(8):: repnl(2),rmunl(2),vavgnl(2)
  complex(8):: smrho(k1,k2,k3,nsp),vxcnl(k1,k2,k3,nsp)
  complex(8):: vxnl(k1,k2,k3,nsp), vcnl(k1,k2,k3,nsp),tpibai

  !      type(s_lat)::  slat
  integer :: ip,i,i1,i2,i3,lcut,ng,np,ipr,n1,n2,n3,ngabc(3),nnn
  real(8),allocatable :: ggr(:,:),agr(:,:),gagr(:,:),rho(:,:)
  real(8),allocatable :: enl(:,:,:),vnl(:,:,:,:),enlbk(:,:,:)
  real(8),allocatable:: dvxcdgr(:,:,:,:),grho2_updn(:,:,:,:),grho(:,:,:,:,:),gv(:,:)
  real(8),allocatable::  grho2_updn_forcall(:,:,:,:)
  complex(8),allocatable:: zgrho(:,:,:),gzgrho(:,:,:),zggrho(:,:)
  complex(8),allocatable:: fgrd(:,:,:),fn(:,:,:),fg(:,:),fgg(:)
  real(8):: alat,vol,xx,fac,tpiba,pi,smmin,sss
  integer ::ig,dummy,j,isp
  logical::  debug=.false., plottest=.false. !newmode=.true. ,
  real(8),allocatable:: r_smrho(:,:,:,:)
  !!== Setup ==
  call tcn('vxcnlm')
  if(debug) print *,'smvxcm:calling vxcnlm sum=',sum(abs(smrho))
  call getpr(ipr)
  ngabc =lat_nabc
  n1= ngabc(1)
  n2= ngabc(2)
  n3= ngabc(3)
  if(abs(n1-k1)+abs(n2-k2)+abs(n3-k3)/=0) call rx('vxcnlm: ni/=ki')
  ng    = lat_ng
  allocate(gv(ng,3))
  call dcopy(3*ng,rv_a_ogv,1,gv,1)
  alat  =lat_alat
  vol   =lat_vol
  np = n1*n2*n3 !k1*k2*k3
  if(debug) print *,'vxcnlm: sum check smrho=',sum(abs(smrho))
  !!== New mode start here ==
  !!== Obtain grho= \nabla smrho (on real mesh) ==
  !  if(newmode) then
     allocate(zgrho(np,3,nsp))
     do  i = 1, nsp
        call grfmsh ( 201 , alat , ng , rv_a_ogv , iv_a_okv , k1 , k2 &
             , k3 , n1 , n2 , n3 , smrho ( 1 , 1 , 1 , i ) , zgrho ( 1 , 1 &
             , i ) , dummy ) !dummy = zzgrho is  not touched for grfmsh(isw=201,...
     enddo
     allocate(grho(k1,k2,k3,3,nsp))
     call dcopy(3*np*nsp,zgrho,2,grho,1) ! grho contains $\nabla$smrho in real space
     deallocate(zgrho)
     !! == grho2_updn = (\nabla smrho) **2 ==
     allocate(grho2_updn(n1,n2,n3,2*nsp-1) )
     do isp=1,nsp
        do i1=1,n1
           do i2=1,n2
              do i3=1,n3
                 grho2_updn(i1,i2,i3,isp) = sum( grho(i1,i2,i3,:,isp)**2 )
                 if(nsp==2) grho2_updn(i1,i2,i3,3) = &
                      sum(grho(i1,i2,i3,1,:))**2 + sum(grho(i1,i2,i3,2,:))**2 + sum(grho(i1,i2,i3,3,:))**2
              enddo
           enddo
        enddo
     enddo
     !!== call xcpbe in abinit ==
     allocate( vnl(k1,k2,k3,nsp) )
     allocate( enl(k1,k2,k3))
     allocate( dvxcdgr(k1,k2,k3,3))
     fac=1d0/2d0  !This fac is required since rp (:,:,isp=1) contains total density in the case of nsp=1.
     if(nsp==2) fac=1d0
     ! i
     allocate(r_smrho(k1,k2,k3,nsp))
     r_smrho=fac*dreal(smrho)
     ! allocate for calling a subroutine
     allocate(grho2_updn_forcall(n1,n2,n3,2*nsp-1))
     do isp=1,2*nsp-1
        do i3=1,n3
           do i2=1,n2
              do i1=1,n1
                 grho2_updn_forcall(i1,i2,i3,isp)=fac**2*grho2_updn(i1,i2,i3,isp)
              enddo
           enddo
        enddo
     enddo
     call xcpbe(exci=enl,npts=n1*n2*n3,nspden=nsp, &
          option=2,&!  Choice of the functional =2:PBE-GGA
          order=1, &!  order=1 means we only calculate first derivative of rho*exc(rho,\nable rho).
                !  rho_updn=fac*dreal(smrho),vxci=vnl,ndvxci=0,ngr2=2*nsp-1,nd2vxci=0,  !Mandatory Arguments
          rho_updn=r_smrho,vxci=vnl,ndvxci=0,ngr2=2*nsp-1,nd2vxci=0,&  ! & Mandatory Arguments
                !  dvxcdgr=dvxcdgr, grho2_updn=fac**2*grho2_updn)   !Optional Arguments
          dvxcdgr=dvxcdgr, grho2_updn=grho2_updn_forcall)   !Optional Arguments
     deallocate(grho2_updn_forcall) ! deallocate temporary data
     deallocate(r_smrho)

     !!=== Output: converted to Ry.===
     enl = 2d0*enl !in Ry.
     vnl = 2d0*vnl !in Ry.
     dvxcdgr= 2d0*dvxcdgr !in Ry.
     !!== vxcnl is given ==
     fac=1d0/2d0
     if(nsp==2) fac=1d0
     allocate(fg(ng,3),fn(k1,k2,k3), fgg(ng) )
     allocate(fgrd(k1,k2,k3))
     do isp=1,nsp
        do j=1,3 !x,y,z components of grho.
           fn(:,:,:) = fac*grho(:,:,:,j,isp)*dvxcdgr(:,:,:,isp)      !exchange part for spin density
           do i1=1,n1
              do i2=1,n2
                 do i3=1,n3
                    fn(i1,i2,i3) = fn(i1,i2,i3) &
                         + sum(grho(i1,i2,i3,j,1:nsp)) * dvxcdgr(i1,i2,i3,3) !correlation part for total density
                    ! sum(grho(:,:,:,j,1:nsp),dim=5) means that a sum only for 1:nsp.-->this cause compile error in gfortran
                    ! n (complex)
                    !          print *,'sumcheck fn j =',sum(abs(fn))
                 enddo
              enddo
           enddo
           call fftz3(fn,n1,n2,n3,k1,k2,k3,1,0,-1)  ! fn (complex) is converted from real space to reciprocal space
           call gvgetf(ng,1,iv_a_okv,k1,k2,k3,fn,fg(1,j))
           !          print *,'sumcheck fg j =',sum(abs(fg(:,j)))
        enddo
        !!== make i G fg ---> FFT back ==
        !         print *, 'check xxx fg(1,1:3) sum =',sum(abs(gv(1,1:3) *fg(1,1:3)))
        pi = 4d0*datan(1d0)
        tpiba = 2d0*pi/alat
        tpibai = dcmplx(0d0,1d0)*tpiba
        do ig=1,ng
           fgg(ig) =  tpibai*sum( gv(ig,1:3) * fg(ig,1:3))
           !           print *, ig,abs(fgg(ig)),sum(gv(ig,1:3)**2),sum(fg(ig,1:3)**2)
        enddo
        !         print *, 'check xxx fgg sum =',sum(abs(fgg)),ng,sum(abs(gv)),sum(abs(fg))
        call gvputf(ng,1,iv_a_okv,k1,k2,k3,fgg,fgrd)
        call fftz3(fgrd,n1,n2,n3,k1,k2,k3,1,0,1)
        vxcnl(:,:,:,isp) = vnl(:,:,:,isp) - dreal(fgrd)
     enddo

     !!=== plottest check write for debug ===
     if(plottest) then
        isp=1
        do i1=1,1
           do i2=1,n2
              do i3=1,n3
                 write(8006,"(3i4,10e12.4)") i1,i2,i3,vxcnl(i1,i2,i3,isp) ,fgrd(i1,i2,i3)
                 write(9006,"(3i4,10e12.4)") i1,i2,i3,enl(i1,i2,i3)
              enddo
              write(8006,*)
              write(9006,*)
           enddo
        enddo
        !        stop ' test end of plot 8006'
     endif
     deallocate(fgrd)
     deallocate(fn,fgg,fg)

     !!=== vxnl and vcnl are dummy now ===
     vxnl=0d0 !dummy now
     vcnl=0d0 !dummy now

     !!== Make reps, rmu ==
     repnl=0d0
     rmunl=0d0
     vavgnl=0d0
     do  i = 1, nsp
        do i3 = 1, n3
           do i2 = 1, n2
              do i1 = 1, n1
                 repnl(i) = repnl(i) + dble(smrho(i1,i2,i3,i))*enl(i1,i2,i3)
                 rmunl(i) = rmunl(i) + dble(smrho(i1,i2,i3,i))*vxcnl(i1,i2,i3,i) !all total
                 vavgnl(i) = vavgnl(i) + vxcnl(i1,i2,i3,i)                       !all total
              enddo
           enddo
        enddo
        repnl(i)  = repnl(i)*vol/(n1*n2*n3)
        rmunl(i)  = rmunl(i)*vol/(n1*n2*n3)
        vavgnl(i) = vavgnl(i)/(n1*n2*n3)
     enddo
     if(plottest) then
        allocate(enlbk(k1,k2,k3))
        enlbk=enl
     endif
     deallocate(grho,grho2_updn,dvxcdgr,vnl,enl)
     call tcx('vxcnlm')
     !!* This is the end of newmode.
     !!--------------------------------------------------------------
end subroutine vxcnlm
   
subroutine smcorm(nbas,sspec,ng,gv,  cgh1,cgh2,lfoc1,lfoc2)
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,ispec
  use m_struc_def           
  use m_lmfinit,only:lat_alat
  use m_lattic,only: lat_vol,rv_a_opos
  !- For foca, add together density of smoothed part of core
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read:
  !i     Stored:    *
  !i     Passed to: corprm
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: alat vol ocy
  !i     Stored:    *
  !i     Passed to: *
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !o Outputs
  !o   cgh1  :Portion of smoothed core that is treated directly
  !o   cgh2  :Portion of smoothed core that is treated perturbatively
  !o   lfoc1 :returned nonzero if any site lfoca is direct (1)
  !o   lfoc2 :returned nonzero if any site lfoca is perturbative
  !u Updates
  !u   02 Jul 05  skip sites for which cofh=0
  !u    1 May 00  Adapted from nfp smc_mesh.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ng,nbas,lfoc1,lfoc2
  real(8):: gv(ng,3)
  type(s_spec)::sspec(*)
  !      type(s_lat)::slat
  double complex cgh1(ng),cgh2(ng)
  ! ... Local parameters
  integer:: k0 , nlmx , kmax , ib , is , lfoc , i
  double precision :: tau(3),v(3),alat,vol,qcorg,qcorh,qsc,cofg,cofh, &
       ceh,rfoc,z
  parameter (k0=3, nlmx=25)
  double complex gkl(0:k0,nlmx)
  alat=lat_alat
  vol=lat_vol
  kmax = 0
  ! --- Accumulate FT of smooth-Hankel foca heads ---
  call dpzero(cgh1,  2*ng)
  call dpzero(cgh2,  2*ng)
  lfoc1 = 0
  lfoc2 = 0
  do  ib = 1, nbas
     is=ispec(ib)
     tau=rv_a_opos(:,ib) 
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     !       qc = qcorg+qcorh
     !        if (iprint() .ge. 50) write(stdo,351) qc,lfoc,qcorg,qcorh
     !  351   format(' qc=',f12.6,'   lfoc',i2,'   qcorg,qcorh',2f12.6)
     if (cofh /= 0) then
        if (lfoc == 1) then
           lfoc1 = 1
           do  i = 1, ng
              v(1) = gv(i,1)
              v(2) = gv(i,2)
              v(3) = gv(i,3)
              call hklft ( v , rfoc , ceh , tau , alat , kmax , 1 , k0 , rv_a_ocy , gkl )
              cgh1(i) = cgh1(i) + cofh*gkl(0,1)/vol
           enddo
        else if (lfoc == 2) then
           call rx('smcorm: we do not allow now lfoc=2 anymore takao') !Aug2010
           lfoc2 = 1
           do  i = 1, ng
              v(1) = gv(i,1)
              v(2) = gv(i,2)
              v(3) = gv(i,3)
              call hklft ( v , rfoc , ceh , tau , alat , kmax , 1 , k0 , rv_a_ocy , gkl )
              cgh2(i) = cgh2(i) + cofh*gkl(0,1)/vol
           enddo
        endif
     endif
  enddo
end subroutine smcorm
end module m_smvxcm
