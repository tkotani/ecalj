!>XC potential for 0th component
module m_smvxcm
  public smvxcm
  private
contains
  subroutine smvxcm(lfrce,smrho,smpot,smvxc,smvx,smvc,smexc,repsm,repsmx,repsmc,rmusm,rvmusm,rvepsm,f) !XC potential for 0th components.
    use m_supot,only: rv_a_ogv,iv_a_okv, lat_ng,n1,n2,n3
    use m_lmfinit,only: rv_a_ocy, lat_alat, lxcfun=>lxcf,nsp,nbas
    use m_lattic,only: lat_vol
    use m_lgunit,only:stdo
    use m_xclda,only: evxcp,evxcv
    use m_ftox
    !i We only allow lxcfun (set at m_lmfinit) are one of
    !         :  1 Ceperly-Alder
    !         :  2  Barth-Hedin (ASW fit)
    !         :  103 PBE-GGA
    !i   lfrce : =1 calculate contribution to forces f
    !i   n1,n2,n3 mesh
    !i   smrho : density   on mesh
    !i   smpot : potential on mesh 
    !i   In addition, n_sH^c contribution is added. See corprm
    !o Outputs
    !o   smvxc :ex-corr  potential of smoothed density + core corrections
    !o   smvx  :exchange potential of smoothed density + core corrections
    !o   smvc  :correlation potential of smoothed density + core corrections
    !o   smexc :ex-corr  energy density of smoothed density + core corrections
    !o   smpot :smooth total potential; smvxc is added to smpot
    !o   repsm :integrated exchange-correlation energy
    !o   repsmx:integrated exchange energy
    !o   repsmc:integrated correlation energy
    !o   rmusm :\int (smrho + n^c_sH,a) * vxc[rhosm+n^c_sH,a]
    !o         :where n^c_sH,a = portion of core treated directly
    !o   rvmusm:int (smrho) * vxc[rhosm+n^c_sH,a]
    !o   rvepsm:int (smrho) * exc[rhosm+n^c_sH,a]
    !o   f     :contribution to forces from xc potential
    !r Remarks
    !   We have to calculate vxc for n0+n^c_sH,a (to take into account the spilout). See Eq.(31) (but typo n_0^Zcv==>n_0^cv).
    ! ----------------------------------------------------------------------
    implicit none
    integer :: lfrce
    real(8):: f(3,nbas) , repsm(nsp) , repsmx(nsp) , repsmc(nsp) , rmusm(nsp) &
         , rvmusm(nsp) , rvepsm(nsp) !, focexc(nsp) , focex(nsp) , focec(nsp) , focvxc(nsp)
    complex(8):: smrho(n1,n2,n3,nsp),smpot(n1,n2,n3,nsp), &
         smvxc(n1,n2,n3,nsp),smvx(n1,n2,n3,nsp),smvc(n1,n2,n3,nsp),smexc(n1,n2,n3)
    integer:: i , n123, ng , lfoc1, lfoc2 , iprint , excsan
    complex(8) ,allocatable :: cgh1_zv(:)
    complex(8) ,allocatable :: dxcv_zv(:)
    real(8) :: vol,sum1,sum2,vxcavg(nsp),x1,x2,alat
    character(180) :: outs
    integer ::iwdummy
    complex(8),allocatable ::smrho_w(:,:,:,:), smcor_w(:),smc(:,:,:)
    integer:: nnn,isp
    real(8):: sss ,srshift,swmin
    real(8),parameter:: minimumrho=1d-14
    logical::enforce_positive_smrho
    logical:: iprx=.true.
    include 'mpif.h'
    integer:: procid=0,ier=0,i1,i2,i3
    integer,parameter::master=0
    call mpi_comm_rank(mpi_comm_world,procid,ier)
    iprx=.false.
    if(procid==master) iprx= .TRUE. 
    call tcn('smvxc')
    ng    = lat_ng
    vol  = lat_vol
    alat = lat_alat
    vol  = lat_vol
    allocate(smc(n1,n2,n3),source=(0d0,0d0))
    allocate(cgh1_zv(ng))
    allocate(smrho_w,mold=smrho)
    call smcorm(nbas,ng,rv_a_ogv,cgh1_zv,lfoc1)
    if(lfoc1 == 1) then
      call gvputf(ng , 1 , iv_a_okv , n1 , n2 , n3 , cgh1_zv , smc ) 
      call fftz3(smc,n1,n2,n3,n1,n2,n3,1,0,1)
    endif
    forall(isp=1:nsp) smrho_w(:,:,:,isp)= smrho(:,:,:,isp)+ smc/nsp !add n^c_sH,a for lfoca=1  smrho_w=n0+n^c_smH,a (for atoms of lfoc1)
    if(enforce_positive_smrho()) then !== negative smrho check== This is also similar with what is done in mkpot.
       nnn   = count(dreal(smrho_w)<minimumrho) 
       if(nnn>0) then ! For GGA, we have to supply positive rho. See enforce_positive_smrho section in mkpot.
          swmin=minval(dreal(smrho_w))
          smrho_w=smrho_w+ minimumrho + abs(swmin)
          if(iprx) write(6,ftox) 'smvxcm: smrho_w<minimumrho  number,min(smrho_w)=',nnn,swmin
       else
          if(iprx) write(6,ftox) 'smvxcm: all smrho_w is positive'
       endif
    endif
    rvmusm = 0d0
    call smvxc2(0, nsp,lxcfun,vol,n1,n2,n3,smrho_w,smvxc,smvx,smvc,smexc,repsm,repsmx,repsmc,rmusm,vxcavg )
    smpot(:,:,:,1:nsp) = smpot(:,:,:,1:nsp)+smvxc(:,:,:,1:nsp)
    do  i = 1, nsp
       rvmusm(i)=dreal(sum(smvxc(:,:,:,i)*smrho(:,:,:,i)))*vol/(n1*n2*n3)
       rvepsm(i)=dreal(sum(smexc(:,:,:)*smrho(:,:,:,i)))*vol/(n1*n2*n3)
    enddo
    deallocate(cgh1_zv)
    if(lfrce /= 0) then !Force from foca sm-head; cgh1 is workspace ---
       allocate(cgh1_zv(ng*nsp))
       f=0d0
       if (lfoc1>0) then
          call fftz3(smvxc,n1,n2,n3,n1,n2,n3,nsp,0,-1)
          call gvgetf(ng, nsp, iv_a_okv,n1,n2,n3,smvxc,cgh1_zv  )
          call smvxc4(nbas, nsp, alat, vol, rv_a_ocy, ng,rv_a_ogv,cgh1_zv,f )
       endif
       if (allocated(cgh1_zv)) deallocate(cgh1_zv)
    endif
    call tcx('smvxc')
  end subroutine smvxcm
  subroutine smvxc2(mode,nsp,lxcfun,vol,n1,n2,n3,smrho, smvxc,smvx,smvc,smexc, rhoeps,rhoex,rhoec,rhomu,vxcavg)
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
    !i         :  1    GGA PBE
    !i   vol   :cell volume
    !i   slat  :struct containing information about the lattice
    !i   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
    !i   n1,n2,n3 dimensions of smrho,smpot for smooth mesh density
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
    integer :: mode,nsp,n1,n2,n3,lxcfun
    real(8) :: rhoeps(nsp),rhoex(nsp),rhoec(nsp),rhomu(nsp), vxcavg(nsp),vol
    complex(8) smvxc(n1,n2,n3,nsp),smrho(n1,n2,n3,nsp), &
         smvx(n1,n2,n3,nsp),smvc(n1,n2,n3,nsp), &
         smexc(n1,n2,n3) !,dsmvxc(n1,n2,n3,nsp)
    ! ... Local parameters
    integer :: i,i1,i2,i3,lxcf,nx,iprint,n1x,lxcg,nglob,mode1
    parameter (n1x=512)
    real(8) :: alfa,dfdr,dvdr,f,f1,f2,rrho,fac,rrmin
    real(8) :: repnl(nsp),rmunl(nsp),vavgnl(nsp)
    real(8) :: vxc2(n1x,2),vxc1(n1x,2), &
         vx1(n1x,2),vc1(n1x,2), &
         exc2(n1x),exc1(n1x), &
         exc1x(n1x),exc1c(n1x)
    real(8) :: rho(n1x),rhos(n1x,2)
    character(180) :: outs
    integer:: ixxx
    real(8):: rhomin
    logical :: newmode=.true.
    integer:: nnn,isp
    real(8):: sss ,smmin(nsp)
    call tcn('smvxc2')
    if (n1 > n1x) call rxi('smvxc2: increase n1x, need',n1)
    lxcf = mod(lxcfun,100)     !1 or 2 or 3
    lxcg = mod(lxcfun/100,100) !1 for GGA
    alfa = 2d0/3d0
    fac = tiny(0d0)**(1d0/3d0) !dmach(1)**(1d0/3d0)
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
    if(lxcg/=0.and.lxcg/=1) call rx('smvxc2: lxcg/=0.and.lxcg/=1')
    if(lxcg == 1) then  !GGA part. Gradient corrections to potential, energy
       call vxcnlm(nsp,n1,n2,n3,smrho, repnl,rmunl,vavgnl,smvxc)
       rhoeps =  repnl
       rhomu  =  rmunl
       vxcavg =  vavgnl
       repnl  = 0d0
       rmunl  = 0d0
       vavgnl = 0d0
    endif
    ! ... LDA, GGA Printout
    if (nx > 0.and.iprint()>=20) write(stdo,ftox)' smvxcm (warning) mesh density ' &
         //'negative at ',nx,'point; rhomin=',rrmin
    if (lxcg /= 0) then
       rhoeps = rhoeps + repnl
       rhomu = rhomu + rmunl
       vxcavg = vxcavg + vavgnl
    endif
    if (iprint() >= 30) then! ... Printout, total potential
       do i=1,nsp
          write(stdo,ftox)'  smvxc2: smooth isp rhoeps rhomu vxcavg=' &
               ,i,ftof(rhoeps(i)),ftof(rhomu(i)),ftof(vxcavg(i))
       enddo
    endif
    call tcx('smvxc2')
  end subroutine smvxc2
  subroutine smvxc3(vol,nsp,n1,n2,n3,smrho,smcor,dsmvxc, smvxc,rmuxcc) !Smooth core density times dvxc/drho
    use m_lgunit,only:stdo
    !i   vol   :cell volume
    !i   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
    !i   n1,n2,n3 dimensions of smrho,smpot for smooth mesh density
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
    integer :: nsp,n1,n2,n3
    real(8) :: rmuxcc(nsp),vol
    complex(8) smvxc(n1,n2,n3,nsp),smcor(n1,n2,n3), &
         dsmvxc(n1,n2,n3,nsp),smrho(n1,n2,n3,nsp)
    integer :: i,i1,i2,i3
    complex(8) cadd,csum(nsp)
    rmuxcc = 0
    csum=0d0
    do  i = 1, nsp
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
  subroutine smvxc4(nbas,nsp,alat,vol,cy,ng,gv,cvxc,f) !- For foca, adds force from shift of smH-head against Vxc.
    use m_lgunit,only:stdo
    use m_struc_def  
    use m_lattic,only: rv_a_opos
    use m_lmfinit,only: ispec
    use m_hansr,only:corprm
    !i   cy    :Normalization constants for spherical harmonics
    !i   ng    :number of G-vectors
    !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
    !i   cvxc  :Fourier transform of smooth vxc potential.
    !o   f     :force from shift of smH-head against Vxc added to f.
    implicit none
    integer :: nbas,nsp,ng
    integer,parameter:: k0=3
    real(8):: gv(ng,3),alat,vol,cy(1),f(3,nbas),fff(3)
    complex(8):: cvxc(ng,nsp),gkl(0:k0,1),ccc,cvxci
    integer :: kmax,ib,is,lfoc,i,kb,iprint
    real(8) :: tau(3),v(3),pi,tpiba,qcorg,qcorh,qsc,cofg, cofh,ceh,rfoc,z,sum1,sum2,sum3,xx
    pi = 4d0*datan(1d0)
    tpiba = 2d0*pi/alat
    kmax = 0
    if (iprint()>= 50) write(stdo,"(/' xc-force from foca:')")
    do  ib = 1, nbas
       is =ispec(ib)
       tau=rv_a_opos(:,ib) 
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
       if (lfoc > 0 .AND. cofh /= 0) then
          fff=0d0
          do  i = 1, ng
             v(:) = gv(i,:)
             call hklft(v,rfoc,ceh,tau,alat,kmax,1,k0,cy,gkl)
             ccc = cofh*gkl(0,1)/vol
             cvxci = 0.5d0 * (cvxc(i,1) + cvxc(i,nsp))
             xx = -dimag(dconjg(cvxci) * ccc)
             fff = fff + xx*gv(i,:)
          enddo
          fff = fff*vol*tpiba
          f(:,ib) = f(:,ib) + fff
          forall(kb = 1:nbas) f(:,kb) = f(:,kb) - fff/nbas
       endif
    enddo
    if(iprint()>=50) write(stdo,"(i4,3f12.6)") (ib,f(1,ib),f(2,ib),f(3,ib),ib = 1,nbas)
  end subroutine smvxc4
  ! Kotani's version newmode with xcpbe.F in abinit Aug2010
  subroutine vxcnlm(nsp,n1,n2,n3,smrho,repnl,rmunl,vavgnl,vxcnl) != Gradient correction to smoothed rho(q) tabulated on a mesh. abinit
    use m_supot,only: iv_a_okv,rv_a_ogv
    use m_xcpbe,  only: xcpbe
    use m_lmfinit,only:    lat_alat
    use m_lattic,only: lat_vol
    use m_supot,only: lat_ng
    ! i   lxcg  : dummy now.  (need to set option in xcpbe)
    ! i   smrho(n1,n2,n3,nsp)
    ! o   repnl : integral smrho * eps
    ! o   rmunl : integral smrho * vxc
    ! o   vavgnl:average NL XC potential
    ! o   vxcnl : XC potential on uniform mesh.
    implicit none
    integer :: lxcg,n1,n2,n3,nsp
    real(8):: repnl(nsp),rmunl(nsp),vavgnl(nsp)
    complex(8):: smrho(n1,n2,n3,nsp),vxcnl(n1,n2,n3,nsp),tpibai
    integer :: ip,i,i1,i2,i3,lcut,ng,np,ipr,nnn
    real(8),allocatable :: ggr(:,:),rho(:,:)
    real(8),allocatable :: enl(:,:,:),vnl(:,:,:,:),enlbk(:,:,:)
    real(8),allocatable:: dvxcdgr(:,:,:,:),grho2_updn(:,:,:,:),grho(:,:,:,:,:),gv(:,:)
    real(8),allocatable::  grho2_updn_forcall(:,:,:,:)
    complex(8),allocatable:: zgrho(:,:,:),gzgrho(:,:,:),zggrho(:,:)
    complex(8),allocatable:: fgrd(:,:,:),fn(:,:,:),fg(:,:),fgg(:)
    real(8):: alat,vol,xx,fac,tpiba,pi,smmin,sss
    integer ::ig,dummy,j,isp
    logical::  debug=.false. !, plottest=.false. 
    real(8),allocatable:: r_smrho(:,:,:,:)
    complex(8):: wdummy(1)
    call tcn('vxcnlm')
    call getpr(ipr)
    if(abs(n1-n1)+abs(n2-n2)+abs(n3-n3)/=0) call rx('vxcnlm: ni/=ki')
    ng    = lat_ng
    allocate(gv(ng,3))
    call dcopy(3*ng,rv_a_ogv,1,gv,1)
    alat  =lat_alat
    vol   =lat_vol
    np = n1*n2*n3 
    if(debug) print *,'vxcnlm: sum check smrho=',sum(abs(smrho))
    !!== Obtain grho= \nabla smrho (on real mesh) ==
    allocate(zgrho(np,3,nsp))
    do  i = 1, nsp
       call grfmsh(201,alat,ng,rv_a_ogv,iv_a_okv,n1,n2,n3, smrho(1,1,1,i),zgrho(1,1,i), wdummy ) !dummy = zzgrho is  not touched for grfmsh(isw=201,...
    enddo
    allocate(grho(n1,n2,n3,3,nsp))
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
    allocate( vnl(n1,n2,n3,nsp) )
    allocate( enl(n1,n2,n3))
    allocate( dvxcdgr(n1,n2,n3,3))
    fac=1d0/2d0!This is required since rp (:,:,isp=1) contains total density in the case of nsp=1.
    if(nsp==2) fac=1d0
    allocate(r_smrho(n1,n2,n3,nsp))
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
    allocate(fg(ng,3),fn(n1,n2,n3), fgg(ng) )
    allocate(fgrd(n1,n2,n3))
    do isp=1,nsp
       do j=1,3 !x,y,z components of grho.
          fn = fac*grho(:,:,:,j,isp)*dvxcdgr(:,:,:,isp)      !exchange part for spin density
          fn=fn+ sum(grho(:,:,:,j,:),dim=4)*dvxcdgr(:,:,:,3) !correlation part for total density
          call fftz3(fn,n1,n2,n3,n1,n2,n3,1,0,-1)  ! from real space to reciprocal space
          call gvgetf(ng,1,iv_a_okv,n1,n2,n3,fn,fg(1,j))  
       enddo
       !!== make i G fg ---> FFT back ==
       pi = 4d0*datan(1d0)
       tpiba = 2d0*pi/alat
       tpibai = dcmplx(0d0,1d0)*tpiba
       fgg =  tpibai*[(sum(gv(ig,1:3)*fg(ig,1:3)),ig=1,ng)]
       call gvputf(ng,1,iv_a_okv,n1,n2,n3,fgg,fgrd)
       call fftz3(fgrd,n1,n2,n3,n1,n2,n3,1,0,1) !fft back to real space
       vxcnl(:,:,:,isp) = vnl(:,:,:,isp) - dreal(fgrd)
    enddo
    deallocate(fgrd)
    deallocate(fn,fgg,fg)
    !!== Make reps, rmu ==
    do  i = 1, nsp
       repnl(i) = sum(dble(smrho(:,:,:,i))*enl(:,:,:))
       rmunl(i) = sum(dble(smrho(:,:,:,i))*vxcnl(:,:,:,i)) !all total
       vavgnl(i)= sum(vxcnl(:,:,:,i))                  !all total
       repnl(i)  = repnl(i)*vol/(n1*n2*n3)
       rmunl(i)  = rmunl(i)*vol/(n1*n2*n3)
       vavgnl(i) = vavgnl(i)/(n1*n2*n3)
    enddo
    deallocate(grho,grho2_updn,dvxcdgr,vnl,enl)
    call tcx('vxcnlm')
  end subroutine vxcnlm
  subroutine smcorm(nbas,ng,gv, cgh1,lfoc1) !For foca, add together density of smoothed part of core
    use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,ispec
    use m_struc_def           
    use m_lmfinit,only:lat_alat
    use m_lattic,only: lat_vol,rv_a_opos
    use m_hansr,only: corprm
    !i   nbas  :size of basis
    !i   ng    :number of G-vectors
    !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
    !o Outputs
    !o   cgh1  :Portion of smoothed core that is treated directly
    !o   lfoc1 :returned nonzero if any site lfoca is direct (1)
    implicit none
    integer,parameter:: k0=3
    integer :: ng,nbas,lfoc1,lfoc2
    real(8):: gv(ng,3)
    complex(8):: cgh1(ng),gkl(0:k0,1)
    integer:: kmax,ib,is,lfoc,i
    real(8) :: tau(3),v(3),alat,vol,qcorg,qcorh,qsc,cofg,cofh, ceh,rfoc,z
    alat=lat_alat
    vol=lat_vol
    kmax = 0
    ! --- Accumulate FT of smooth-Hankel foca heads ---
    cgh1=0d0
    lfoc1 = 0
    do  ib = 1, nbas
       is=ispec(ib)
       tau=rv_a_opos(:,ib) 
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
       if (cofh /= 0.and. lfoc == 1) then
          lfoc1 = 1
          do  i = 1, ng
             v = gv(i,:)
             call hklft ( v,rfoc,ceh,tau,alat,kmax,1,k0,rv_a_ocy,gkl )
             cgh1(i) = cgh1(i) + cofh*gkl(0,1)/vol
          enddo
       endif
    enddo
  end subroutine smcorm
  subroutine grfmsh(isw,alat,ng,gv,kv,n1,n2,n3,fn, fgrd,flap) !- Gradient and Laplacian of a function tabulated on a uniform mesh
  !i Inputs
  !i   isw   :1s digit
  !i         : 0 input function fn is real
  !i         : 1 input function fn is complex
  !i         :10s digit
  !i         : 0 Input fn is in real space
  !i         : 1 Input fn is in reciprocal space
  !i         :100s digit handles output
  !i         : 0 No output
  !i         : 1 Output fgrd, in reciprocal space, no flap
  !i         : 2 Output fgrd, in real space, no flap
  !i         : 3 Output flap, in reciprocal space, no fgrd
  !i         : 4 Output flap, in real space, no fgrd
  !i         : 5 Output fgrd and flap, in reciprocal space
  !i         : 6 Output fgrd and flap, in real space
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   vol   :cell volume
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !i   kv    :indices for gather/scatter operations (gvlist.f)
  !i   n1..n3:dimensions of fn,fgrd,flap
  !i   n1..n3:size of uniform mesh for which fn,fgrd,flap are tabulated
  !i   fn    :function on uniform mesh, either in real or recip space,
  !i         :depending on isw
  !o Outputs
  !o   fgrd  :gradient of fn either in real or recip space, depending on isw
  !o   flap  :laplacian of fn either in real or recip space, depending on isw
  !l Local variables
  !l   gi    :i * G (for gradient)
  !l   g2    :(i * G)^2 (for Laplacian)
  !l   fg    :gradient of function, G space
  !l   fl    :Laplacian of function, G space
  !r Remarks
  !r Sample test:
  !r   cmfft -f9f15.10 -248 rho | mc -f9f15.10  . igx -xe | cmfft -f9f15.10 -i -248 .
  !u Updates
  !u   08 Apr 09 First created
  ! ----------------------------------------------------------------------
 implicit none
 integer :: ng,isw,n1,n2,n3
 integer :: kv(ng,3)
 real(8) :: alat,gv(ng,3)
 complex(8) fn(n1,n2,n3),fgrd(n1,n2,n3,3),flap(n1,n2,n3)
 integer :: i,isw0,isw1,isw2
 real(8) :: pi,tpiba,g2
 complex(8):: tpibai,gi(3)
 complex(8),allocatable:: fg(:),fgg(:,:),fg2(:)
 call tcn('grfmsh')
 pi   = 4d0*datan(1d0)
 tpiba=2*pi/alat
 isw0 = mod(isw,10)
 isw1 = mod(isw/10,10)
 isw2 = mod(isw/100,10)
 if (isw0 /= 1) call rx('grfmsh: gradient of real function not implemented')
 ! ... FT of smooth function to reciprocal space
 if (isw1 == 0) call fftz3(fn,n1,n2,n3,n1,n2,n3,1,0,-1)
 ! ... Gather function G coefficients
 allocate(fg(ng))
 call gvgetf(ng,1,kv,n1,n2,n3,fn,fg)
 ! ... Restore given function to real space
 if (isw1 == 0) call fftz3(fn,n1,n2,n3,n1,n2,n3,1,0,1)
 ! ... Make iG * f(G)  and (iG)^2 * f(G)
 fg(1) = 0
 tpibai = dcmplx(0d0,1d0)*tpiba
 allocate(fgg(ng,3),fg2(ng))
 do  i = 1, ng
    gi(:) = tpibai*(gv(i,:))
    g2 = -tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
    fgg(i,:) = gi(:) * fg(i)
    fg2(i)   = g2    * fg(i)
 enddo
 ! ... Scatter gradient, laplacian into respective arrays
 if (isw2 == 1 .OR. isw2 == 2 .OR. isw2 == 5 .OR. isw2 == 6) then
    call gvputf(ng,1,kv,n1,n2,n3,fgg(1,1),fgrd(1,1,1,1))
    call gvputf(ng,1,kv,n1,n2,n3,fgg(1,2),fgrd(1,1,1,2))
    call gvputf(ng,1,kv,n1,n2,n3,fgg(1,3),fgrd(1,1,1,3))
 endif
 if (isw2 == 3 .OR. isw2 == 4 .OR. isw2 == 5 .OR. isw2 == 6) then
    call gvputf(ng,1,kv,n1,n2,n3,fg2,     flap)
 endif
 ! ... Gradient, laplacian in real space
 if (isw2 == 2 .OR. isw2 == 6) then
    call fftz3(fgrd(1,1,1,1),n1,n2,n3,n1,n2,n3,1,0,1)
    call fftz3(fgrd(1,1,1,2),n1,n2,n3,n1,n2,n3,1,0,1)
    call fftz3(fgrd(1,1,1,3),n1,n2,n3,n1,n2,n3,1,0,1)
 endif
 if (isw2 == 4 .OR. isw2 == 6) call fftz3(flap,n1,n2,n3,n1,n2,n3,1,0,1)
 deallocate(fg,fgg,fg2)
 call tcx('grfmsh')
end subroutine grfmsh
end module m_smvxcm
