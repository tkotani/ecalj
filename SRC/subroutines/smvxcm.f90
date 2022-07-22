module m_smvxcm
  public smvxcm
contains
  subroutine smvxcm(lfrce,smrho,smpot,smvxc,smvx,smvc,smexc,repsm,repsmx,repsmc,rmusm,rvmusm,rvepsm,f)
    use m_supot,only: rv_a_ogv,iv_a_okv, lat_nabc, lat_ng,k1,k2,k3
    use m_lmfinit,only: rv_a_ocy,   lat_alat, lxcf,nsp,nbas
    use m_lattic,only: lat_vol
    use m_lgunit,only:stdo
    use m_xclda,only: evxcp,evxcv
    !- XC potential for smooth mesh density
    ! ----------------------------------------------------------------------
    !i Inputs
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
    !o   f     :contribution to forces from xc potential
    !l Local variables
    !l  lxcfun :1s digit sets local xc functional
    !l         :  1    Ceperly-Alder
    !l         :  2    Barth-Hedin (ASW fit)
    !l         :  103 PBE-GGA
    !r Remarks
    !r   smoothed core is partition into core1 + core2.  All atoms with
    !r   lfoc1=1 belong to core1
    !r  *core1 is included directly into smrho; the nonlinear XC potential
    !r   is computed from vxc[smrho+smcor1].
    !rxxx  *core2 removed!
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
    integer :: lfrce,ngabc(3),lxcfun
    equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
    real(8):: f(3,nbas) , repsm(nsp) , repsmx(nsp) , repsmc(nsp) , rmusm(nsp) &
         , rvmusm(nsp) , rvepsm(nsp) !, focexc(nsp) , focex(nsp) , focec(nsp) , focvxc(nsp)
    complex(8):: smrho(k1,k2,k3,nsp),smpot(k1,k2,k3,nsp), &
         smvxc(k1,k2,k3,nsp),smvx(k1,k2,k3,nsp),smvc(k1,k2,k3,nsp),smexc(k1,k2,k3)
    integer:: i , k123 , n1 , n2 , n3 , ng , lfoc1, lfoc2 , iprint , excsan
    complex(8) ,allocatable :: cgh1_zv(:)
    complex(8) ,allocatable :: cgh2_zv(:)
    complex(8) ,allocatable :: dxcv_zv(:)
    double precision :: vol,sum1,sum2,vxcavg(nsp),x1,x2,alat
    character(180) :: outs
    integer ::iwdummy
    complex(8),allocatable ::smrho_w(:), smcor_w(:)
    integer:: nnn,isp
    real(8):: sss ,srshift,swmin
    real(8),parameter:: minimumrho=1d-14
    logical::enforce_positive_smrho
    logical:: iprx=.true.
    include 'mpif.h'
    integer:: procid=0,ier=0
    integer,parameter::master=0
    call mpi_comm_rank(mpi_comm_world,procid,ier)
    iprx=.false.
    if(procid==master) iprx= .TRUE. 
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
    call smcorm(nbas,ng,rv_a_ogv,cgh1_zv,cgh2_zv, lfoc1,lfoc2)
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
          if(sss<minimumrho) then 
             nnn=nnn+1
             if(sss<swmin) then
                swmin=sss
             endif
          endif
       enddo
       if(nnn>0) then
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
    rvmusm = 0d0
    ! --- Direct branch (lfoc2 .eq. 0) ---
    if (lfoc2 == 0) then
       call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho_w,sum1,sum2)
       call smvxc2 ( 0 ,  nsp , lxcfun , vol , n1 , n2 , n3 ,&
            k1 , k2 , k3 , smrho_w , smvxc , smvx , smvc , smexc  &
            , repsm , repsmx , repsmc , rmusm , vxcavg )
       smpot(:,:,:,1:nsp) = smpot(:,:,:,1:nsp)+smvxc(:,:,:,1:nsp)
       do  i = 1, nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),smrho(1,1,1,i),rvmusm(i),x2)
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),  smrho(1,1,1,i),rvepsm(i),x2)
       enddo
    endif
    if (lfoc2 /= 0) call rx('smvxcm: not suppor perturbation branch now')
    deallocate(cgh2_zv,cgh1_zv)
    ! --- Force from foca sm-head; cgh1 is workspace ---
    if (lfrce /= 0) then
       allocate(cgh1_zv(ng*nsp))
       call dpzero(f,3*nbas)
       if (lfoc1 > 0 .OR. lfoc2 > 0) then
          call fftz3(smvxc,n1,n2,n3,k1,k2,k3,nsp,0,-1)
          call gvgetf(ng, nsp, iv_a_okv , k1 , k2 , k3 , smvxc , cgh1_zv  )
          call smvxc4(nbas, nsp, alat, vol, rv_a_ocy, ng , rv_a_ogv , cgh1_zv , f )
       endif
       if (allocated(cgh1_zv)) deallocate(cgh1_zv)
    endif
    call tcx('smvxc')
  end subroutine smvxcm
  subroutine smvxc2(mode,nsp,lxcfun,vol,n1,n2,n3,k1,k2,k3,smrho, &
       smvxc,smvx,smvc,smexc,       rhoeps,rhoex,rhoec,rhomu,vxcavg)
    use m_ftox
    use m_lgunit,only:stdo
    use m_xclda,only: evxcv,evxcp
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
    double precision :: rhoeps(nsp),rhoex(nsp),rhoec(nsp),rhomu(nsp), vxcavg(nsp),vol
    double complex smvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp), &
         smvx(k1,k2,k3,nsp),smvc(k1,k2,k3,nsp), &
         smexc(k1,k2,k3) !,dsmvxc(k1,k2,k3,nsp)
    ! ... Local parameters
    integer :: i,i1,i2,i3,lxcf,nx,iprint,n1x,lxcg,nglob,mode1
    parameter (n1x=512)
    double precision :: alfa,dfdr,dvdr,f,f1,f2,rrho,fac,dmach,rrmin
    double precision :: repnl(nsp),rmunl(nsp),vavgnl(nsp)
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
    real(8):: sss ,smmin(nsp)
    call tcn('smvxc2')
    if (n1 > n1x) call rxi('smvxc2: increase n1x, need',n1)
    lxcf = mod(lxcfun,100)
    lxcg = mod(lxcfun/100,100)
    alfa = 2d0/3d0
    fac = dmach(1)**(1d0/3d0)
    rhoeps = 0
    rhoex  = 0
    rhoec  = 0
    rhomu  = 0
    vxcavg = 0
    rrmin = 0
    nx = 0
    mode1 = mod(mode/10,10)
    smvxc=0d0
    smvx=0d0
    smvc=0d0
    smexc=0d0
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
             call dcopy(n1,vxc1(1,i),1,smvxc(1,i2,i3,i),2)
             call dcopy(n1,vx1(1,i),1,smvx(1,i2,i3,i),2)
             call dcopy(n1,vc1(1,i),1,smvc(1,i2,i3,i),2)
          enddo
          call dcopy(n1,exc1,1,smexc(1,i2,i3),2)
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
    rhoeps = rhoeps*f
    rhoex = rhoex*f
    rhoec = rhoec*f
    rhomu = rhomu*f
    vxcavg = vxcavg/(n1*n2*n3)
    ! ... Gradient corrections to potential, energy
    if(lxcg/=0.and.lxcg/=1) call rx('smvxc2: lxcg/=0.and.lxcg/=1')
    if (lxcg == 1) then  !        print *,'smvxcm:calling vxcnlm sum smrho=',sum(abs(smrho))
       call vxcnlm(lxcg,nsp,k1,k2,k3,smrho, repnl,rmunl,vavgnl,smvx,smvc,smvxc)
       rhoeps =  repnl
       rhomu  =  rmunl
       vxcavg =  vavgnl
       repnl  = 0d0
       rmunl  = 0d0
       vavgnl = 0d0
    endif
    ! ... LDA, GGA Printout
    if (nx > 0) call info5(20,0,0,' smvxcm (warning) mesh density ' &
         //'negative at %i point%?#n>1#s##:  rhomin=%;3g',nx,nx,rrmin,0,0)
    if (lxcg /= 0) then
       rhoeps = rhoeps + repnl
       rhomu = rhomu + rmunl
       vxcavg = vxcavg + vavgnl
    endif
    if (iprint() >= 30) then! ... Printout, total potential
       do i=1,nsp
          write(stdo,ftox)'smooth isp rhoeps rhomu vxcavg=' &
               ,i,ftof(rhoeps(i)),ftof(rhomu(i)),ftof(vxcavg(i))
       enddo
    endif
    call tcx('smvxc2')
  end subroutine smvxc2
  subroutine smvxc3(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,smcor,dsmvxc, smvxc,rmuxcc)
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
    integer :: nsp,k1,k2,k3,n1,n2,n3
    double precision :: rmuxcc(nsp),vol
    double complex smvxc(k1,k2,k3,nsp),smcor(k1,k2,k3), &
         dsmvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp)
    integer :: i,i1,i2,i3
    double complex cadd,csum(nsp)
    rmuxcc = 0
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
  end subroutine smvxc3
  subroutine smvxc4(nbas,nsp,alat,vol,cy,ng,gv,cvxc,f)
    use m_lgunit,only:stdo
    use m_struc_def  
    use m_lattic,only: rv_a_opos
    use m_lmfinit,only: ispec
    !- For foca, adds force from shift of smH-head against Vxc.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   nsp   :number of spin channels
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
    integer :: nbas,nsp,ng
    real(8):: gv(ng,3) , alat , vol , cy(1) , f(3,nbas)
    double complex cvxc(ng,nsp)
    integer :: k0,nlmx,kmax,ib,is,lfoc,i,kb,iprint
    double precision :: tau(3),v(3),pi,tpiba,qcorg,qcorh,qsc,cofg, &
         cofh,ceh,rfoc,z,sum1,sum2,sum3,xx
    parameter (k0=3, nlmx = 9)
    double complex gkl(0:k0,nlmx),ccc,cvxci
    pi = 4d0*datan(1d0)
    tpiba = 2d0*pi/alat
    kmax = 0
    ! --- Loop over sites ---
    if (iprint() >= 50) write(stdo,400)
    do  ib = 1, nbas
       is=ispec(ib)
       tau=rv_a_opos(:,ib) 
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
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
    if (iprint() >= 50) write(stdo,340) (ib,f(1,ib),f(2,ib),f(3,ib),ib = 1,nbas)
340 format(i4,3f12.6)
400 format(/' xc-force from foca:')
  end subroutine smvxc4
  !!= Gradient correction to smoothed rho(q) tabulated on a mesh =
  !!*Kotani's version newmode with xcpbe.F in abinit Aug2010
  subroutine vxcnlm(lxcg,nsp,k1,k2,k3,smrho,repnl,rmunl,vavgnl,vxnl,vcnl,vxcnl)
    use m_supot,only: iv_a_okv,rv_a_ogv
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
    real(8):: repnl(nsp),rmunl(nsp),vavgnl(nsp)
    complex(8):: smrho(k1,k2,k3,nsp),vxcnl(k1,k2,k3,nsp)
    complex(8):: vxnl(k1,k2,k3,nsp), vcnl(k1,k2,k3,nsp),tpibai
    integer :: ip,i,i1,i2,i3,lcut,ng,np,ipr,n1,n2,n3,ngabc(3),nnn
    real(8),allocatable :: ggr(:,:),agr(:,:),gagr(:,:),rho(:,:)
    real(8),allocatable :: enl(:,:,:),vnl(:,:,:,:),enlbk(:,:,:)
    real(8),allocatable:: dvxcdgr(:,:,:,:),grho2_updn(:,:,:,:),grho(:,:,:,:,:),gv(:,:)
    real(8),allocatable::  grho2_updn_forcall(:,:,:,:)
    complex(8),allocatable:: zgrho(:,:,:),gzgrho(:,:,:),zggrho(:,:)
    complex(8),allocatable:: fgrd(:,:,:),fn(:,:,:),fg(:,:),fgg(:)
    real(8):: alat,vol,xx,fac,tpiba,pi,smmin,sss
    integer ::ig,dummy,j,isp
    logical::  debug=.false. !, plottest=.false. !newmode=.true. ,
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
    fac=1d0/2d0!This is required since rp (:,:,isp=1) contains total density in the case of nsp=1.
    if(nsp==2) fac=1d0
    allocate(r_smrho(k1,k2,k3,nsp))
    r_smrho=fac*dreal(smrho)
    ! allocate for calling a subroutine
    allocate(grho2_updn_forcall(n1,n2,n3,2*nsp-1))
    do isp=1,2*nsp-1
       grho2_updn_forcall(:,:,:,isp)=fac**2*grho2_updn(:,:,:,isp)
    enddo
    call xcpbe(exci=enl,npts=n1*n2*n3,nspden=nsp, &
         option=2,&!  Choice of the functional =2:PBE-GGA
         order=1, &!  order=1 means we only calculate first derivative of rho*exc(rho,\nable rho).
         rho_updn=r_smrho,vxci=vnl,ndvxci=0,ngr2=2*nsp-1,nd2vxci=0,&  ! & Mandatory Arguments
         dvxcdgr=dvxcdgr, grho2_updn=grho2_updn_forcall)   !Optional Arguments
    deallocate(grho2_updn_forcall) 
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
          fn = fac*grho(:,:,:,j,isp)*dvxcdgr(:,:,:,isp)      !exchange part for spin density
          fn=fn+ sum(grho(:,:,:,j,:),dim=4)*dvxcdgr(:,:,:,3) !correlation part for total density
          call fftz3(fn,n1,n2,n3,k1,k2,k3,1,0,-1)  ! from real space to reciprocal space
          call gvgetf(ng,1,iv_a_okv,k1,k2,k3,fn,fg(1,j))  
       enddo
       !!== make i G fg ---> FFT back ==
       pi = 4d0*datan(1d0)
       tpiba = 2d0*pi/alat
       tpibai = dcmplx(0d0,1d0)*tpiba
       fgg =  tpibai*[(sum(gv(ig,1:3)*fg(ig,1:3)),ig=1,ng)]
       call gvputf(ng,1,iv_a_okv,k1,k2,k3,fgg,fgrd)
       call fftz3(fgrd,n1,n2,n3,k1,k2,k3,1,0,1) !fft back to real space
       vxcnl(:,:,:,isp) = vnl(:,:,:,isp) - dreal(fgrd)
    enddo
    !!=== plottest check write for debug ===
    ! if(plottest) then
    !    isp=1
    !    do i1=1,1
    !       do i2=1,n2
    !          do i3=1,n3
    !             write(8006,"(3i4,10e12.4)") i1,i2,i3,vxcnl(i1,i2,i3,isp) ,fgrd(i1,i2,i3)
    !             write(9006,"(3i4,10e12.4)") i1,i2,i3,enl(i1,i2,i3)
    !          enddo
    !          write(8006,*)
    !          write(9006,*)
    !       enddo
    !    enddo
    ! endif
    deallocate(fgrd)
    deallocate(fn,fgg,fg)
    !!=== vxnl and vcnl are dummy now ===
    vxnl=0d0 !dummy now
    vcnl=0d0 !dummy now
    !!== Make reps, rmu ==
    do  i = 1, nsp
       repnl(i) = sum(dble(smrho(:,:,:,i))*enl(:,:,:))
       rmunl(i) = sum(dble(smrho(:,:,:,i))*vxcnl(:,:,:,i)) !all total
       vavgnl(i)= sum(vxcnl(:,:,:,i))                  !all total
       repnl(i)  = repnl(i)*vol/(n1*n2*n3)
       rmunl(i)  = rmunl(i)*vol/(n1*n2*n3)
       vavgnl(i) = vavgnl(i)/(n1*n2*n3)
    enddo
    !if(plottest) then
    !   allocate(enlbk(k1,k2,k3))
    !   enlbk=enl
    !endif
    deallocate(grho,grho2_updn,dvxcdgr,vnl,enl)
    call tcx('vxcnlm')
  end subroutine vxcnlm
  subroutine smcorm(nbas,ng,gv,  cgh1,cgh2,lfoc1,lfoc2)
    use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,ispec
    use m_struc_def           
    use m_lmfinit,only:lat_alat
    use m_lattic,only: lat_vol,rv_a_opos
    !- For foca, add together density of smoothed part of core
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
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
    integer :: ng,nbas,lfoc1,lfoc2
    real(8):: gv(ng,3)
    double complex cgh1(ng),cgh2(ng)
    integer:: k0 , nlmx , kmax , ib , is , lfoc , i
    double precision :: tau(3),v(3),alat,vol,qcorg,qcorh,qsc,cofg,cofh, ceh,rfoc,z
    parameter (k0=3, nlmx=25)
    double complex gkl(0:k0,nlmx)
    alat=lat_alat
    vol=lat_vol
    kmax = 0
    ! --- Accumulate FT of smooth-Hankel foca heads ---
    cgh1=0d0
    cgh2=0d0
    lfoc1 = 0
    lfoc2 = 0
    do  ib = 1, nbas
       is=ispec(ib)
       tau=rv_a_opos(:,ib) 
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
       if (cofh /= 0) then
          if (lfoc == 1) then
             lfoc1 = 1
             do  i = 1, ng
                v = gv(i,:)
                call hklft ( v , rfoc , ceh , tau , alat , kmax , 1 , k0 , rv_a_ocy , gkl )
                cgh1(i) = cgh1(i) + cofh*gkl(0,1)/vol
             enddo
          else if (lfoc == 2) then
             call rx('smcorm: we do not allow now lfoc=2 anymore takao') !Aug2010
             lfoc2 = 1
             do  i = 1, ng
                v = gv(i,:)
                call hklft ( v , rfoc , ceh , tau , alat , kmax , 1 , k0 , rv_a_ocy , gkl )
                cgh2(i) = cgh2(i) + cofh*gkl(0,1)/vol
             enddo
          endif
       endif
    enddo
  end subroutine smcorm
end module m_smvxcm
