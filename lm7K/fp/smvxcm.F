Cchp1 cmdvar= real(8),rv_p_,ogv
Cchp1 cmdvar= real(8),rv_p_,ocy
      subroutine smvxcm(ssite,sspec,slat,nbas,lfrce,k1,k2,k3,smrho,
     .smpot,smvxc,smvx,smvc,smexc,repsm,repsmx,repsmc,rmusm,rvmusm,
     .rvepsm,focexc,focex,focec,focvxc,f)

      use m_struc_def  !Cgetarg

C- XC potential for smooth mesh density
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct containing site-specific information
Ci     Elts read: spec pos
Ci     Stored:
Ci     Passed to: smcorm smvxc4
Ci   sspec :struct containing species-specific information
Ci     Elts read:
Ci     Stored:
Ci     Passed to: smcorm smvxc4 corprm
Ci   slat  :struct containing information about the lattice
Ci     Elts read: nabc ng ogv okv vol alat ocy
Ci     Stored:
Ci     Passed to: smcorm
Ci   nbas  :size of basis
Ci   lfrce :1, calculate contribution to forces
Ci   k1,k2,k3 dimensions smooth crystal densities, potentials on a mesh
Ci   smrho :smooth valence density on uniform mesh
Co Outputs
Co   smvxc :ex-corr  potential of smoothed density + core corrections
Co   smvx  :exchange potential of smoothed density + core corrections
Co   smvc  :correlation potential of smoothed density + core corrections
Co   smexc :ex-corr  energy density of smoothed density + core corrections
Co   smpot :smooth total potential; smvxc is added to smpot
Co   repsm :integrated exchange-correlation energy
Co   repsmx:integrated exchange energy
Co   repsmc:integrated correlation energy
Co   rmusm :int (smrho + smcor1) * vxc[rhosm+smcor1]
Co         :where smcor1 = portion of core treated directly
Co   rvmusm:int (smrho) * vxc[rhosm+smcor1]
Co   rvepsm:int (smrho) * exc[rhosm+smcor1]
Co   focexc:FOCA exchange-correlation energy:
Co         :int (smcor2) * vxc[rhosm+smcor1]
Co         :where smcor2 = portion of core treated perturbatively
Co   focex :exchange part of focexc
Co   focec :correlation part of focexc
Co   focvxc:integral of FOCA exchange-correlation potential:
Co         :int (smcor2) * (smrho) * (dvxc/drho)
Co   f     :contribution to forces from xc potential
Cl Local variables
Cl  lxcfun :1s digit sets local xc functional
Cl         :  1    Ceperly-Alder
Cl         :  2    Barth-Hedin (ASW fit)
Cl         :  3,4  LD part of PW91 and PBE
Cl         :100s digit sets gradient corrections
Cl         :  0    LSDA
Cl         :  1    Langreth-Mehl
Cl         :  2    PW91
Cl         :  3    PBE
Cl         :  4    PBE with Becke exchange
Cr Remarks
Cr   smoothed core is partition into core1 + core2.  All atoms with
Cr   lfoc1=1 belong to core1; those with lfoc2=1 belong to core2.
Cr  *core1 is included directly into smrho; the nonlinear XC potential
Cr   is computed from vxc[smrho+smcor1].
Cr  *core2 is included perturbatively: its contribution to the vxc
Cr   is computed from the expansion
Cr     vxc[rho + smcor2] = vxc[rho] + smcor2 * (dvxc/drho)
Cr                       = vxc[rho] + dvxc
Cr   The perturbation correction to int (smrho * vxc) is then
Cr     focvxc = int smrho * smcor2 * (dvxc/drho)
Cr   If the perturbation approach is exact,
Cr     (focvxc+rvmusm) -> rvmusm when computed with smcor2=0
Cr   The corresponding XC energy density is
Cr     exc[rho + smcor2] = exc[rho] + smcor2 * (dexc/drho)
Cr                       = exc[rho] + smcor2 * (vxc-exc)/rho
Cr   The perturbation correction to the XC energy is then
Cr     int smcor2 * (vxc-exc) = focexc - int smcor2 exc[smrho+smcor1]
Cu Updates
Cu   21 Apr 09 Handles GGA functionals
Cu   02 Jul 05 skip sites for which cofh=0
Cu   25 Jun 04 return smexc,rvepsm
Cu   14 Jun 02 rhoex and rhoec (T. Miyake)
Cu    8 Feb 02 smvx and smvc (T. Miyake)
Cu   13 Jun 00 spin polarized
Cu    1 May 00 Adapted from nfp sxc_smooth.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,lfrce,k1,k2,k3,ngabc(3)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      real(8):: f(3,nbas) , repsm(2) , repsmx(2) , repsmc(2) , rmusm(2) 
     ., rvmusm(2) , rvepsm(2) , focexc(2) , focex(2) , focec(2) , focvxc(2)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex smrho(k1,k2,k3,2),smpot(k1,k2,k3,2),
     .smvxc(k1,k2,k3,2),smvx(k1,k2,k3,2),smvc(k1,k2,k3,2),
     .smexc(k1,k2,k3)
C ... Local parameters
Cchp1       integer:: stdo , i , k123 , n1 , n2 , n3 , lgunit , ng , ogv 
Cchp1      ., okv , ocy , lfoc1 , lfoc2 , iprint , excsan 
Cchp1      ., lxcfun , nglob , nsp
Cchp1        integer:: stdo , i , k123 , n1 , n2 , n3 , lgunit , ng , ogv 
Cchp1      . , okv , lfoc1 , lfoc2 , iprint , excsan , lxcfun , nglob , nsp 
        integer:: stdo , i , k123 , n1 , n2 , n3 , lgunit , ng , okv 
     .  , lfoc1 , lfoc2 , iprint , excsan , lxcfun , nglob , nsp 
        real(8),pointer :: rv_p_ogv(:)

       real(8),pointer :: rv_p_ocy(:)

      complex(8) ,allocatable :: cgh1_zv(:)
      complex(8) ,allocatable :: cgh2_zv(:)
      complex(8) ,allocatable :: dxcv_zv(:)

      double precision vol,sum1,sum2,vxcavg(2),x1,x2,alat
      character*180 outs
C ... Heap
      integer w(1)
      common /w/ w
      complex(8),allocatable ::smrho_w(:), smcor_w(:) 

      stdo = lgunit(1)
      nsp  = nglob('nsp')
      lxcfun = nglob('lxcf')
      call tcn('smvxc')

c      i_copy_size=size(slat%nabc)
c      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      ngabc=slat%nabc
      ng=slat%ng
Cchp1       ogv=slat%ogv
       rv_p_ogv => slat%rv_p_ogv 

      okv=slat%okv
      vol=slat%vol
      alat=slat%alat
      vol=slat%vol
Cchp1       ocy=slat%ocy
       rv_p_ocy => slat%rv_p_ocy 


c      if (excsan(lxcfun,0) .lt. 0) x1 = excsan(lxcfun,stdo)

C ... Sum of foca hankel heads; break into direct and pert. parts
c      call defcc (osmrho,-k1*k2*k3*nsp)
      allocate(smrho_w(k1*k2*k3*nsp))
      smrho_w=0d0
      allocate(cgh1_zv(ng), cgh2_zv(ng))
Cchp1       call smcorm ( nbas , ssite , sspec , slat , ng , w ( ogv ) , 
Cchp1      .cgh1_zv , cgh2_zv , lfoc1 , lfoc2 )
       call smcorm ( nbas , ssite , sspec , slat , ng , rv_p_ogv , cgh1_zv 
     . , cgh2_zv , lfoc1 , lfoc2 ) 


C ... smrho_w = smrho + smoothed core from foca hankel heads
      k123 = 2*k1*k2*k3
      if (lfoc1 .eq. 1) then
        call gvputf ( ng , 1 , w ( okv ) , k1 , k2 , k3 , cgh1_zv 
     .  , smrho_w ) !w ( osmrho ) )
        call fftz3(smrho_w,n1,n2,n3,k1,k2,k3,1,0,1)
        if (nsp .eq. 2) then
          call dscal(k123,.5d0,smrho_w,1)
          call dpscop(smrho_w,smrho_w,k123,1,1+k123,1d0)
        endif
        call dpadd(smrho_w,smrho,1,k123*nsp,1d0)
      else
        call dpcopy(smrho,smrho_w,1,k123*nsp,1d0)
      endif

C ... Force density strictly positive definite
C      print *, 'density strictly pos def?'
C      call smrpos(w(osmrho),k1,k2,k3,n1,n2,n3)
      rvmusm(2) = 0d0
      focexc(2) = 0d0
      focex(2)  = 0d0
      focec(2)  = 0d0
      focvxc(2) = 0d0

C --- Direct branch (lfoc2 .eq. 0) ---
      if (lfoc2 .eq. 0) then
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho_w,sum1,sum2)
        call smvxc2(0,slat, nsp,lxcfun,vol,n1,n2,n3,k1,k2,k3,
     .  smrho_w,smvxc,smvx,smvc,smexc,w,
     .  repsm,repsmx,repsmc,rmusm,vxcavg)
        call dpadd(smpot,smvxc,1,2*k1*k2*k3*nsp,1d0)
        do  i = 1, nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),
     .    smrho(1,1,1,i),rvmusm(i),x2)
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),
     .    smrho(1,1,1,i),rvepsm(i),x2)
        enddo
        focexc(1) = 0d0
        focex(1)  = 0d0
        focec(1)  = 0d0
        focvxc(1) = 0d0
      endif

C --- Perturbation branch (lfoc2 .ne. 0) ---
      if (lfoc2 .ne. 0) then
        allocate(dxcv_zv(k1*k2*k3*nsp))
        call smvxc2 ( 1 , slat,nsp , lxcfun , vol , n1 , n2 , n3 , k1 , k2 
     .  , k3 , smrho_w , smvxc , smvx , smvc , smexc , dxcv_zv 
     .  , repsm , repsmx , repsmc , rmusm , vxcavg )
        call dpadd(smpot,smvxc,1,2*k1*k2*k3*nsp,1d0)
        do  i = 1, nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),
     .    smrho(1,1,1,i),rvmusm(i),x2)
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),
     .    smrho(1,1,1,i),rvepsm(i),x2)
        enddo
        deallocate(smrho_w)
C   ... Assemble core tails for linearized treatment on mesh
C       w(osmcor) = portion of core treated perturbatively
c        osmcor = osmrho
        allocate(smcor_w(k1*k2*k3*nsp))
        call gvputf ( ng , 1 , w ( okv ) , k1 , k2 , k3 , cgh2_zv 
     .  , smcor_w) !w ( osmcor ) )

        call fftz3(smcor_w,n1,n2,n3,k1,k2,k3,1,0,1)
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smcor_w,sum1,sum2)
        call dpzero(focexc,2)
        call dpzero(focex,2)
        call dpzero(focec,2)
        do  i = 1, nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvx(1,1,1,i),smcor_w,
     .    x1,x2)
          focex(i)  = x1/nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvc(1,1,1,i),smcor_w,
     .    x1,x2)
          focec(i)  = x1/nsp
          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smvxc(1,1,1,i),smcor_w,
     .    x1,x2)
          focexc(i) = x1/nsp
C         Add this term to focexc to make focexc=pert corr to rhov*exc
C          call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smexc(1,1,1),w(osmcor),
C     .      x1,x2)
        enddo
C       Peturbation correction to smvxc
        call smvxc3 ( vol , nsp , n1 , n2 , n3 , k1 , k2 , k3 , smrho 
     .  , smcor_w , dxcv_zv , smvxc , focvxc )

        call dpadd(smpot,smvxc,1,2*k1*k2*k3*nsp,1d0)
        if (iprint() .ge. 30) then
          outs = ' '
          call awrit8('%x   foca'//
     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N   foca#%2j#'//
     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',
     .    outs,len(outs),0,
     .    focexc(1)+focexc(2),nsp,focexc,focexc(2),
     .    focvxc(1)+focvxc(2),nsp,focvxc,focvxc(2))
          call awrit1('%a  charge  =%;12,6D',outs,len(outs),-stdo,sum1)
        endif
        deallocate(smcor_w,dxcv_zv)
      endif
      deallocate(cgh2_zv,cgh1_zv)
c      call rlse(osmrho)

C --- Force from foca sm-head; cgh1 is workspace ---
      if (lfrce .ne. 0) then
        allocate(cgh1_zv(ng*nsp))
        call dpzero(f,3*nbas)
        if (lfoc1.gt.0 .or. lfoc2.gt.0) then
          call fftz3(smvxc,n1,n2,n3,k1,k2,k3,nsp,0,-1)
          call gvgetf ( ng , nsp , w ( okv ) , k1 , k2 , k3 , smvxc , cgh1_zv 
     .    )

Cchp1           call smvxc4 ( nbas , nsp , ssite , sspec , alat , vol , w ( ocy 
Cchp1      .    ) , ng , w ( ogv ) , cgh1_zv , f )
Cchp1            call smvxc4 ( nbas , nsp , ssite , sspec , alat , vol , rv_p_ocy 
Cchp1      .     , ng , w ( ogv ) , cgh1_zv , f ) 
            call smvxc4 ( nbas , nsp , ssite , sspec , alat , vol , rv_p_ocy 
     .      , ng , rv_p_ogv , cgh1_zv , f ) 



        endif
        if (allocated(cgh1_zv)) deallocate(cgh1_zv)

      endif

      call tcx('smvxc')

      end subroutine smvxcm


      subroutine smvxc2(mode,slat,nsp,lxcfun,vol,n1,n2,n3,k1,k2,k3,smrho,
     .smvxc,smvx,smvc,smexc,dsmvxc,rhoeps,rhoex,rhoec,rhomu,vxcavg)
      use m_struc_def, only: s_lat
C- Makes smooth part of xc potential smvxc and optionally dsmvxc/drho
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do not make dsmvxc/drho
Ci         :1 make dsmvxc/drho
Ci         :10s digit
Ci         : 0 calculated LDA for density as is.
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci   nsp   :number of spin channels
Ci  lxcfun :switch defining xc functional (evxcv.f)
Ci         :1s digit sets local xc functional
Ci         :  1    Ceperly-Alder
Ci         :  2    Barth-Hedin (ASW fit)
Ci         :  3,4  LD part of PW91 and PBE
Ci         :100s digit sets gradient corrections
Ci         :  0    LSDA
Ci         :  1    Langreth-Mehl
Ci         :  2    PW91
Ci         :  3    PBE
Ci         :  4    PBE with Becke exchange
Ci   vol   :cell volume
Ci   slat  :struct containing information about the lattice
Ci   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   smrho :smooth density on uniform mesh
Co Outputs
Co   smvxc :xc potential of smoothed density (no core contr.)
Co   smvx  :exchange potential of smoothed density + core corrections
Co   smvc  :correlation potential of smoothed density + core corrections
Co   dsmvxc:dvxc/drho (mode=1)
Co   rhoeps:integrated exchange-correlation energy
Co   rhoex :integrated exchange energy
Co   rhoec :integrated correlation energy
Co   rhomu :integrated exchange-correlation potential
Co   vxcavg:average xc potential
Cr Remarks
Cr   For perturbation treatment, take numerical derivatives
Cr   df/dr = d/dr (vxc*r**alfa) instead of d/dr vxc because
Cr   f is nearly linear for alpha=2/3.
Cr
Cr   In the spin polarized case, the smooth core density is not
Cr   spin polarized.  Thus to calc. vxc(rho+drho, m+dm) - vxc(rho,m)
Cr   we use dm=0 and thus drho1 = drho2 = drho/2; thus
Cr     dvxc = lim_drho->0  vxc(rho+drho,rho1+drho/2) - vxc(rho,rho1)
Cr
Cu Updates
Cu   20 Nov 09 New 10s digit for mode
Cu   21 Apr 09 Handles GGA functionals
Cu   14 Jun 02 rhoex and rhoec (T. Miyake)
Cu    8 Feb 02 smvx and smvc (T. Miyake)
Cu   12 Jun 00 spin polarized
Cu    1 May 00 Adapted from nfp vxcd_smooth.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nsp,k1,k2,k3,n1,n2,n3,lxcfun
      double precision rhoeps(2),rhoex(2),rhoec(2),rhomu(2),
     .vxcavg(2),vol
      double complex smvxc(k1,k2,k3,2),smrho(k1,k2,k3,2),
     .smvx(k1,k2,k3,2),smvc(k1,k2,k3,2),
     .smexc(k1,k2,k3,2),dsmvxc(k1,k2,k3,2)
C ... Local parameters
      integer i,i1,i2,i3,lgunit,lxcf,nx,stdo,iprint,n1x,lxcg,nglob,mode1
      parameter (n1x=512)
      double precision alfa,dfdr,dvdr,f,f1,f2,rrho,fac,dmach,rrmin
      double precision repnl(2),rmunl(2),vavgnl(2)
      double precision vxc2(n1x,2),vxc1(n1x,2),
     .vx1(n1x,2),vc1(n1x,2),
     .exc2(n1x),exc1(n1x),
     .exc1x(n1x),exc1c(n1x)
      double precision rho(n1x),rhos(n1x,2)
      character*180 outs

      integer(4):: ixxx
      real(8):: rhomin
      type(s_lat)::slat
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      smrho(:,:,:,2)=0d0
c      write(6,"(a,12d13.5)")'rrrr: smchk smrho 1'
c     &  , sum(smrho(:,:,:,1)), sum(abs(smrho(:,:,:,1)))
c      write(6,"(a,12d13.5)")'rrrr: smchk smrho 2'
c     &  , sum(smrho(:,:,:,2)), sum(abs(smrho(:,:,:,2)))
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,'smvxc2: xxx'

      if (n1 .gt. n1x) call rxi('smvxc2: increase n1x, need',n1)
      stdo = lgunit(1)
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

      call dpzero(smvxc,2*k1*k2*k3*nsp)
      call dpzero(smvx,2*k1*k2*k3*nsp)
      call dpzero(smvc,2*k1*k2*k3*nsp)
      call dpzero(smexc,2*k1*k2*k3)

C     Vector of points for each i2,i3.
      do  i3 = 1, n3
        do  i2 = 1, n2

C     ... rho = total rho for this vec; rhos = spin pol rho
          call dcopy(n1,smrho(1,i2,i3,1),2,rho,1)
          call dcopy(n1,smrho(1,i2,i3,1),2,rhos,1)
          if (nsp .eq. 2) then
            call dcopy(n1,smrho(1,i2,i3,2),2,rhos(1,2),1)
            call daxpy(n1,1d0,rhos(1,2),1,rho,1)
          endif

C     ... Put df/dr into dsmvxc
          if (mod(mode,10) .ne. 0) then
            if (lxcg .ne. 0) call rx('smvxc3: perturbation treatment '//
     .        'not implemented with GGA')
C           Add rho*fac/2 into rhos+, rhos- and fac*rho into rho
            do  i = 1, nsp
              call daxpy(n1,fac/2,rho,1,rhos(1,i),1)
            enddo
            call dscal(n1,1+fac,rho,1)
C           Exchange potential at rho+drho
            do  i = 1, nsp
              call evxcv(rho,rhos(1,i),n1,1,lxcf,
     .        exc2,exc1x,exc1c,vxc2(1,i),vx1(1,i),vc1(1,i))
            enddo
            if (mode1 .eq. 2) then
              do  i1 = 1, n1
                if (rhos(i1,1).lt.0 .or. rhos(i1,nsp).lt.0) then
                  vxc2(i1,1) = 0
                  vxc2(i1,nsp) = 0
                endif
              enddo
            endif
C           Restore rho,rhos; also add -drho to rho and -drho/2 to rhos
            call dscal(n1,1/(1+fac),rho,1)
            do  i = 1, nsp
              call daxpy(n1,-fac,rho,1,rhos(1,i),1)
            enddo
            call dscal(n1,(1-fac),rho,1)
C           Exchange potential at rho-drho
            do  i = 1, nsp
              call evxcv(rho,rhos(1,i),n1,1,lxcf,
     .        exc1,exc1x,exc1c,vxc1(1,i),vx1(1,i),vc1(1,i))
            enddo
            if (mode1 .eq. 2) then
              do  i1 = 1, n1
                if (rhos(i1,1).lt.0 .or. rhos(i1,nsp).lt.0 .or.
     .              vxc2(i1,1).eq.0 .or. vxc2(i1,nsp) .eq. 0) then
                  vxc1(i1,1) = 0
                  vxc1(i1,nsp) = 0
                  vxc2(i1,1) = 0
                  vxc2(i1,nsp) = 0
                endif
              enddo
            endif

C           Restore rho,rhos
            call dscal(n1,1/(1-fac),rho,1)
            do  i = 1, nsp
              call daxpy(n1,fac/2,rho,1,rhos(1,i),1)
            enddo

            do  i = 1, nsp
              do  i1 = 1, n1
                if (rho(i1) .gt. 0) then
                  f1 = vxc1(i1,i)*(rho(i1)*(1-fac))**alfa
                  f2 = vxc2(i1,i)*(rho(i1)*(1+fac))**alfa
                  dfdr = (f2-f1)/(2d0*fac*rho(i1))
                  vxc2(i1,i) = dfdr
                else
                  vxc2(i1,i) = 0
                endif
              enddo
            enddo
          endif

C     ... Exchange into smvxc
          if (lxcf .eq. 3 .or. lxcf .eq. 4) then
            call evxcp(rhos(1,1),rhos(1,2),n1,nsp,lxcf,exc1x,exc1c,exc1,
     .        vx1(1,1),vx1(1,2),vc1(1,1),vc1(1,2),vxc1(1,1),vxc1(1,2))
          else
            do  i = 1, nsp
              call evxcv(rho,rhos(1,i),n1,nsp,lxcf,exc1,exc1x,exc1c,
     .          vxc1(1,i),vx1(1,i),vc1(1,i))
            enddo
          endif
          if (mode1 .eq. 2) then
            do  i = 1, nsp
            do  i1 = 1, n1
              if (rho(i1).le.0 .or. rhos(i1,i).le.0) then
                vxc1(1,i) = 0
              endif
            enddo
          enddo
        endif

          do  i = 1, nsp
c            call evxcv(rho,rhos(1,i),n1,nsp,lxcf,exc1,
c     .      exc1x,exc1c,vxc1(1,i),vx1(1,i),vc1(1,i))
            call dcopy(n1,vxc1(1,i),1,smvxc(1,i2,i3,i),2)
            call dcopy(n1,vx1(1,i),1,smvx(1,i2,i3,i),2)
            call dcopy(n1,vc1(1,i),1,smvc(1,i2,i3,i),2)
            call dcopy(n1,exc1,1,smexc(1,i2,i3,1),2)
          enddo

C     ... Perturbation dv/dr into dsmvxc
          if (mod(mode,10) .ne. 0) then
            do  i = 1, nsp
              do  i1 = 1, n1
                rrho = rho(i1)
                if (rrho .gt. 0) then
                  f = vxc1(i1,i) * rrho**alfa
                  dvdr = (vxc2(i1,i) - alfa*f/rrho) / rrho**alfa
                  dsmvxc(i1,i2,i3,i) = dvdr
                else
                  dsmvxc(i1,i2,i3,i) = 0
                endif
              enddo
            enddo
          endif

C     ... Add to integrals
          do  i = 1, nsp
            do  i1 = 1, n1
              rrho = rhos(i1,i)
              rrmin = min(rrho,rrmin)
              if (rrho .lt. 0d0) nx = nx+1
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

C      do  i = 1, nsp
C        call zprm3('LDA smvxc (isp=%i)',i,smvxc(1,1,1,i),n1,n2,n3)
C        call zprm3('dsmvxc (isp=%i)',i,dsmvxc(1,1,1,i),n1,n2,n3)
C      enddo
c      print *,' goto vxcnlm lxcg=',lxcg,sum(abs(smrho)),sum(smrho)

C ... Gradient corrections to potential, energy
      if (lxcg .ne. 0) then
        call vxcnlm(lxcg,nsp,k1,k2,k3,slat,smrho,repnl,rmunl,vavgnl,
     .    smvx,smvc,smvxc)
C        repnl = 0 ; rmunl = 0 ; vavgnl = 0; print *, '!!'
      endif
C      do  i = 1, nsp
C        call zprm3('LDA+GGA smvxc (isp=%i)',i,smvxc(1,1,1,i),n1,n2,n3)
C      enddo

C     call setpr(30)

C ... LDA, GGA Printout
      if (nx .gt. 0) call info5(20,0,0,' smvxcm (warning) mesh density '
     . //'negative at %i point%?#n>1#s##:  rhomin=%;3g',nx,nx,rrmin,0,0)

      if (iprint() .ge. 30 .and. lxcg .ne. 0) then
        call awrit8('%x sm LDA'//
     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%1f#%2j#sm LDA'//
     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
     .    0,rhoeps(1)+rhoeps(2),nsp,rhoeps,rhoeps(2),rhomu(1)+rhomu(2),
     .    nsp,rhomu,rhomu(2))
        call awrit5('%a%?#n==2#%N%3f#  #avg LDA vxc ='//
     .    '%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
     .    -stdo,nsp,(vxcavg(1)+vxcavg(nsp))/2,nsp,vxcavg,vxcavg(2))

        call awrit8('%x sm GGA'//
     .    ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%1f#%2j#sm GGA'//
     .    '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
     .    0,repnl(1)+repnl(2),nsp,repnl,repnl(2),rmunl(1)+rmunl(2),
     .    nsp,rmunl,rmunl(2))
        call awrit5('%a%?#n==2#%N%3f#  #avg GGA vxc ='//
     .    '%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
     .    -stdo,nsp,(vavgnl(1)+vavgnl(nsp))/2,nsp,vavgnl,vavgnl(2))
      endif

      if (lxcg .ne. 0) then
        do  i = 1, nsp
          rhoeps(i) = rhoeps(i) + repnl(i)
C         rhoex(i) = rhoex(i) +
C         rhoec(i) = rhoec(i) +
          rhomu(i) = rhomu(i) + rmunl(i)
          vxcavg(i) = vxcavg(i) + vavgnl(i)
        enddo
      endif

C ... Printout, total potential
      if (iprint() .ge. 30) then
        call awrit8('%x smooth'//
     .  ' rhoeps =%;12,6D %?#n==2#(%;11,6D,%;11,6D)%N%7f#%2j#'//
     .  '  rhomu =%;12,6D %?#n==2#(%;11,6D,%;11,6D)#%2j#',outs,120,
     .  0,rhoeps(1)+rhoeps(2),nsp,rhoeps,rhoeps(2),rhomu(1)+rhomu(2),
     .  nsp,rhomu,rhomu(2))
        call awrit5('%a%?#n==2#%N%7f#  #'//
     .  'avg vxc =%;12,6D %?#n==2#(%;11,6D,%;11,6D)',outs,len(outs),
     .  -stdo,nsp,(vxcavg(1)+vxcavg(nsp))/2,nsp,vxcavg,vxcavg(2))
      endif

C      call zprm3('smvxc',0,smvxc,n1,n2,n3)
C      call zprm3('dsmvxc',0,dsmvxc,n1,n2,n3)
C      if (nsp .eq. 2) then
C        call zprm3('smvxc spin 2',0,smvxc(1,1,1,2),n1,n2,n3)
C        call zprm3('dsmvxc spin 2',0,dsmvxc(1,1,1,2),n1,n2,n3)
C      endif

      end subroutine smvxc2


      subroutine smvxc3(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,smcor,dsmvxc,
     .smvxc,rmuxcc)

C- Smooth core density times dvxc/drho
C ----------------------------------------------------------------------
Ci Inputs
Ci   vol   :cell volume
Ci   n1,n2,n3 uniform mesh on which smrho,smcor,cmvxc defined
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   smrho :smooth density on n1,n2,n3 mesh
Ci   smcor :smooth core density on n1,n2,n3 mesh
Ci   dsmvxc:dvxc/drho on n1,n2,n3 mesh mesh
Co Outputs
Co   smvxc :(dvxc/drho * smcor) = pert. correction to expansion
Co         : vxc[rho + rhoc] = vxc[rho] + rhoc * dvxc/drho
Co   rmuxcc:integral smrho * (dvxc/drho * smcor)
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,k1,k2,k3,n1,n2,n3
      double precision rmuxcc(nsp),vol
      double complex smvxc(k1,k2,k3,nsp),smcor(k1,k2,k3),
     .dsmvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp)
C ... Local parameters
      integer i,i1,i2,i3
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

C     write(stdo,862) csum
C 862 format(' csum=',2f14.8)

      end subroutine smvxc3

      subroutine smvxc4(nbas,nsp,ssite,sspec,alat,vol,cy,ng,gv,cvxc,f)

      use m_struc_def  !Cgetarg

C- For foca, adds force from shift of smH-head against Vxc.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: corprm
Ci   cy    :Normalization constants for spherical harmonics
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   cvxc  :Fourier transform of smooth vxc potential.
Co Outputs
Co   f     :force from shift of smH-head against Vxc added to f.
Cr Remarks
Cu Updates
Cu   02 Jul 05  skip sites for which cofh=0
Cu    1 May 00  Adapted from nfp smc_force.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,ng
      real(8):: gv(ng,3) , alat , vol , cy(1) , f(3,nbas)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      double complex cvxc(ng,nsp)
C ... Local parameters
      integer stdo,k0,nlmx,lgunit,kmax,ib,is,lfoc,i,kb,iprint
      double precision tau(3),v(3),pi,tpiba,qcorg,qcorh,qsc,cofg,
     .cofh,ceh,rfoc,z,sum1,sum2,sum3,xx
      parameter (k0=3, nlmx = 9)
      double complex gkl(0:k0,nlmx),ccc,cvxci

      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      tpiba = 2d0*pi/alat
      kmax = 0

C --- Loop over sites ---
      if (iprint() .ge. 50) write(stdo,400)
      do  ib = 1, nbas

        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pos)
        call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1)

        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        if (lfoc .gt. 0 .and. cofh .ne. 0) then
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
      if (iprint() .ge. 50)
     .write(stdo,340) (ib,f(1,ib),f(2,ib),f(3,ib),ib = 1,nbas)
  340 format(i4,3f12.6)
  400 format(/' xc-force from foca:')

      end subroutine smvxc4

