Cdwdef1 have_common_w_in mkpot , #_of_w_access= 0
Cchp3 cmdvar= type(s_rv1),sv_p_,orhoat
Cchp3 cmdvar= type(s_rv1),sv_p_,oppi
Cchp3 cmdvar= type(s_rv1),sv_p_,otau
Cchp3 cmdvar= type(s_rv1),sv_p_,osig
Cchp3       subroutine mkpot(nbas,ssite,sspec,slat,sham,spot,lfrce,lcplxp,
Cchp3      .k1,k2,k3,smrho,orhoat,qbg,smpot,qmom,vconst,vesrmt,osig,otau,
Cchp3      .oppi,ppn,hab,vab,sab,qval,qsc,gpot0,vval,fes,job,
Cchp3      .vorb,nlibu,lmaxu,lldau,  itest)
Cchp3        subroutine mkpot ( nbas , ssite , sspec , slat , sham , spot 
Cchp3      . , lfrce , lcplxp , k1 , k2 , k3 , smrho , orhoat , qbg , smpot 
Cchp3      . , qmom , vconst , vesrmt , sv_p_osig , otau , oppi , ppn , hab 
Cchp3      . , vab , sab , qval , qsc , gpot0 , vval , fes , job , vorb , 
Cchp3      . nlibu , lmaxu , lldau , itest ) 
Cchp3         subroutine mkpot ( nbas , ssite , sspec , slat , sham , spot 
Cchp3      .  , lfrce , lcplxp , k1 , k2 , k3 , smrho , orhoat , qbg , smpot 
Cchp3      .  , qmom , vconst , vesrmt , sv_p_osig , sv_p_otau , oppi , ppn 
Cchp3      .  , hab , vab , sab , qval , qsc , gpot0 , vval , fes , job , vorb 
Cchp3      .  , nlibu , lmaxu , lldau , itest ) 
Cchp3          subroutine mkpot ( nbas , ssite , sspec , slat , sham , spot 
Cchp3      .   , lfrce , lcplxp , k1 , k2 , k3 , smrho , orhoat , qbg , smpot 
Cchp3      .   , qmom , vconst , vesrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
Cchp3      .   , ppn , hab , vab , sab , qval , qsc , gpot0 , vval , fes , job 
Cchp3      .   , vorb , nlibu , lmaxu , lldau , itest ) 
          subroutine mkpot ( nbas , ssite , sspec , slat , sham , spot 
     .    , lfrce , lcplxp , k1 , k2 , k3 , smrho , sv_p_orhoat , qbg , 
     .    smpot , qmom , vconst , vesrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
     .    , ppn , hab , vab , sab , qval , qsc , gpot0 , vval , fes , job 
     .    , vorb , nlibu , lmaxu , lldau , itest ) 





      use m_struc_def  !Cgetarg

C- Make the potential from the density.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci     Passed to: rhomom smves smvxcm locpot
Ci   sspec :struct containing species-specific information
Ci     Elts read:
Ci     Passed to: rhomom smves smvxcm locpot
Ci   slat  :struct containing information about the lattice
Ci     Elts read: nabc vol
Ci     Passed to: smves smvxcm locpot
Ci   sham  :struct containing information about the hamiltonian
Ci     Elts read: eterms
Ci     Stored:    eterms
Ci   spot  :struct containing information about the potential
Ci     (not used at present)
Ci   lfrce :nonzero =>  contribution to forces
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   k1,k2,k3 dimensions of smrho for smooth crystal density
Ci   smrho :smooth crystal density, on a uniform mesh
Ci   orhoat:local atomic densities (true and smooth parts)
Ci   qbg   :homogeneous background charge
Ci   job   :1s digit
Ci         : 0 stops after energy terms
Ci         : 1 makes potpars also
Ci         :10s digit
Ci         :1 suppress updating potential used to make potpars
Ci         :100s digit
Ci         :1 exclude exchange-correlation potential
Ci         :1000s digit
Ci         :1 Make rveps and rvvxc
Ci         :10000s digit
Ci         :1 write sphere density for each site ib to file rhoMT.ib
Ci ... The following are LDA+U inputs
Ci   vorb  :orbital dependent potential
Ci   nlibu :number of U blocks  used to dimension vorb
Ci   lmaxu :max l for U blocks  used to dimension vorb
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat beginning at dmats(*,lldau(ib))
Co Outputs:
Co   smpot :smooth potential on a uniform mesh:
Co         :Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
Co   qmom  :multipole moments of valence sphere densities
Co   vconst:additional constant potential term
Co   vesrmt:electrostatic potential at rmt
Co   qval  :total valence charge, including semicore states
Co   qsc   :total charge from semicore states (local orbitals)
Co   osig,otau,oppis: augmentation matrices
Co   ppn   :nmto-like potential parameters
Co   hab,vab,sab: augmentation matrices for local parts of the ham.
Co         :See Remarks in augmat.f for their generation and doc.
Co   vval  :coffs to YL expansion of es potential at MT boundary
Co   gpot0 :integrals of gaussians times electrostatic potential
Co   fes   :contribution to the force from electrostatic + xc potential
Co   sham->eterms various integrals for the total energy are stored:
Co         :(1)  ehar   --- not touched here
Co         :(2)  eks    --- not touched here
Co         :(3)  utot   = total electrostatic energy
Co         :(4)  valves = valence rho * estat potential
Co         :(5)  cpnves = core+nuc * estat potential
Co         :(6)  rhoexc = rho * exc
Co         :(7)  rhovxc = rho * vxc
Co         :(8)  sumec  --- not touched here
Co         :(9)  sumtc  --- not touched here
Co         :(10) xcore  = rhoc * total potential
Co         :(11) valvef = smrhov * vsm + sum_ib valvef_ib
Co                        valvef_ib = rhov * vtrue - smrhov * vsm)_ib
Co         :(12) sumt0  --- not touched here
Co         :(13) dq1    --- not touched here
Co         :(14) dq2    --- not touched here
Co         :(15) amom   = system magnetic moment
Co         :(16) sumev  = sphere sum-of-eigenvalues --- not touched here
Co         :(17) rinvxt --- not touched here
Co         :(18) rouvxt --- not touched here
Co         :(19) rhosig --- not touched here
Cl Local variables
Cl   rvmusm:int rhosm * vxc(rhosm+smcor1) where smcor1 is portion
Cl         :of smooth core density treated nonperturbatively.
Cl   rvepsm:int rhosm * exc(rhosm+smcor1) where smcor1 is portion
Cl         :of smooth core density treated nonperturbatively.
Cl   rvepsv:integral of valence density times exc(valence density)
Cl   rvvxcv:integral of valence density times vxc(valence density)
Cl   rhvsm :integral n0~ phi0~
Cl   sgp0  :compensating gaussians * sm-Ves = int (n0~-n0) phi0~
Cl   valfsm:rhvsm - sgp0 + rvmusm + fcvxc0 = n0 (phi0~ + Vxc(n0))
Cl         :NB sgp0 associated with the local parts, as is done in
Cl         :making the ppi matrix elements.
Cl   valftr:local contribution to density * veff
Cl         := sgp0 - fcvxca + sum_ib valvfa_ib
Cl   valvfa:generated by locpot, called valvef there.  Sum_ib of:
Cl         :vefv1 - vefv2
Cl         := int[rho1*(v1-2*Z/r) - (rho2*v2)] - (n0~-n0)*gpotb - focvxc
Cl   lso   :nonzero => spin orbit coupling
Cr Remarks
Cr *The total density is a sum of three terms,
Cr
Cr    n0(mesh) + sum_RL (n_RL(r) - n0_RL(r))
Cr
Cr  The first term is the smooth density on a mesh of points; the
Cr  second is the true density and is defined on a radial mesh for each
Cr  sphere; the last is the 1-center expansion of the smooth density on
Cr  the radial mesh.  (Note: because of l-truncation, n0_R(r) is not
Cr  identical to the one-center expansion of n0(mesh).  The sum of the
Cr  three terms converges rapidly with l because errors in n_R(r) are
Cr  mostly canceled by errors in n0_R(r).)
Cr
Cr *Computation of the electrostatic energy:
Cr  We add and subtract a set of compensating gaussian orbitals
Cr
Cr    n0 + sum_RL Q_RL g_RL + sum_RL (n_RL(r) - n0_RL(r) - Q_RL g_RL)
Cr
Cr  which render the integral of the local part (the last 3 terms)
Cr  zero in each RL channel.  The g_RL must be localized enough that
Cr  their spillout beyond the MT radius is negligible.
Cr
Cr  We define
Cr
Cr    n0~ = n0 + compensating gaussians
Cr
Cr  In the interstitial, the electrostatic potential of n0~ is the true
Cr  estat potential.  The potential of n0 is called phi0 and the
Cr  potential of n0~ is called phi0~.  The total electrostatic energy
Cr  is computed as
Cr    the electrostatic energy of  n0~  +
Cr    the electrostatic energy of (neutral) local parts
Cr
Cr  The first term is computed in subroutine smves;
Cr  the second term is computed in subroutine locpot.
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   27 Apr 05 Added LDA+U stuff (Lambrecht)
Cu   24 Dec 04 Changes for full L.S coupling
Cu    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
Cu   15 Jul 04 (Chantis) Matrix elements Lz.Sz for spin-orbit coupling
Cu   14 Jan 02 rvexv and rvecv (T. Miyake)
Cu   17 Sep 01 Returns qsc.  Altered argument list.
Cu   24 Aug 01 Extended to local orbitals.  Altered argument list.
Cu             Local potentials now have correct boundary conditions
Cu   15 Aug 01 Generates rvepsv and rvvxcv.  Changed call to locpot.
Cu   20 Apr 01 Generates vesrmt
Cu   18 Apr 01 Added ability to exclude exchange-correlation potential
Cu   20 Feb 01 Added ppn to potential parameters generated
Cu   15 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp mk_potential.f
C ------------------------------------------------------------
C to do:
C 1. check that rhov*vxc approximately the same in locpot as smvxc
C 2. ditto for rhov*exc; need finish making rhov*exc in smvxc
C 3. ? See about understanding and changing foxexc both in locpot
C    and smvxc.  The total energy seems to come out just right,
C    but focexc isn't real correction to rhov*exc.  Does it matter?
C    Maybe can't use lfoca=2 is all.
C 4. enable 1000s digit job to pass through
C 5. change locpot so that rvvxcv etc is make only if appropriate
C    1000s job bit set.
C 6. Make <vxc_nn> and compare to integrals made here.
C    If problems, Maybe even look at sm part only, or local part only.
C     implicit none
C ... Passed parameters
      integer nbas,job,lfrce,lcplxp,k1,k2,k3,n0,nppn,nab,i1,i2,i3
      integer nlibu,lmaxu,lldau(nbas)
      parameter (n0=10,nppn=12,nab=9)
      double precision qbg,qval,qsc,vconst,vesrmt(nbas)
Cchp3       integer osig(nbas),otau(nbas),oppi(nbas),orhoat(3,nbas)
Cchp3        integer:: otau(nbas) , oppi(nbas) , orhoat(3,nbas) 
Cchp3         integer:: oppi(nbas) , orhoat(3,nbas) 
Cchp3          integer:: orhoat(3,nbas) 
          type(s_rv1) :: sv_p_orhoat(3,nbas)

         type(s_rv1) :: sv_p_oppi(nbas)

        type(s_rv1) :: sv_p_otau(nbas)

       type(s_rv1) :: sv_p_osig(nbas)

      real(8):: qmom(1) , fes(3,1) , gpot0(1) , vval(1) , ppn(nppn,n0,nbas) 
     ., hab(nab,n0,nbas) , sab(nab,n0,nbas) , vab(nab,n0,nbas)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_ham)::sham
      type(s_pot)::spot

      double complex smrho(k1,k2,k3,*),smpot(k1,k2,k3,*)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,nlibu)
C ... Local parameters
      logical cmdopt,lgors
      character*80 outs
      integer i,lgunit,nglob,ipl,ipr,iprint,n1,n2,n3,ngabc(3),stdo,
     .stdl,lxcfun,lso,isw,nsp,isum
      real(8) ,allocatable :: fxc_rv(:)
      real(8) ,allocatable :: hpot0_rv(:)
      complex(8) ,allocatable :: smvxc_zv(:)
      complex(8) ,allocatable :: smvx_zv(:)
      complex(8) ,allocatable :: smvc_zv(:)
      complex(8) ,allocatable :: smexc_zv(:)

      double precision amom,cpnvsa,cpnves,dq,eh,eks,focexc,
     .focvxc,qsmc,smq,smag,sum2,rhoexc,rhoex,rhoec,rhovxc,rhvsm,sgp0,
     .sqloc,sqlocc,saloc,sumec,sumtc,uat,usm,utot,valfsm,valftr,
     .rvepva,rvexva,rvecva,rvvxva,rvepsa,rvvxca,valvfa,vvesat,valvef,
     .valves,vol,vsum,xcore,zsum,zvnsm,rvepsv(2),rvexv(2),rvecv(2),
     .rvvxcv(2),fcexc0(2),rveps(2),rvvxc(2),fcex0(2),fcec0(2),
     .fcexca(2),fcexa(2),fceca(2),rvmusm(2),rvepsm(2),rmusm(2),
     .vxcavg(2),fcvxc0(2),fcvxca(2),repat(2),repatx(2),repatc(2),
     .rmuat(2),repsm(2),repsmx(2),repsmc(2),rhobg,pi
C     double precision tmp(2)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))


C     This ordering must match sham->eterms; see uham
      double precision eterms(20)
      equivalence (eterms(1),eh)
      equivalence (eterms(2),eks)
      equivalence (eterms(3),utot)
      equivalence (eterms(4),valves)
      equivalence (eterms(5),cpnves)
      equivalence (eterms(6),rhoexc)
      equivalence (eterms(7),rhovxc)
      equivalence (eterms(8),sumec)
      equivalence (eterms(9),sumtc)
      equivalence (eterms(10),xcore)
      equivalence (eterms(11),valvef)
      equivalence (eterms(15),amom)

C ... Heap
Cdwdef1       integer w(1)
Cdwdef1       common /w/ w

ccccccccccccccccccccc
      integer:: itest,wdummy
      complex(8),allocatable:: smpotbk(:,:,:),smpotbkx(:,:,:)
cccccccccccccccccccccc


C      print *, '!!'
C      job = job + 30000

      call tcn('mkpot')
ccccccccccccccccccccccccccccccccc
c      print *,'vvv mkpot sum smrho=',sum(smrho(:,:,:,1)),sum(abs(smrho(:,:,:,1)))
c      open(1898,file='test1vvv')
c      nnnew=0
c      do  i3 = 1, k3
c      do  i2 = 1, k2
c      do  i1 = 1, k1
c        if (dreal(smrho(i1,i2,i3,1)) .lt. 0) nnnew = nnnew+1
cc        write(1898,"(3i5,2d13.6)")i1,i2,i3,dreal(smrho(i1,i2,i3,1))
c      enddo
c      enddo
c      enddo
cc      close(1898)
c      print *,'vvv mkpot sum smrho=',sum(smrho(:,:,:,1)),sum(abs(smrho(:,:,:,1))),nnnew
cc      stop 'xxxxxxxxxx 1111'
ccccccccccccccccccccccccccccccccc

      ipr = iprint()
      ipl = ipr
      lso = isw ( iand(4,int(sham%lncol)) .ne.0 ) + 2 * isw ( iand(32,int(sham%lncol)) .ne.0 
     .) + 3 * isw ( iand(64,int(sham%lncol)) .ne.0 )

      if (lso .ne. 0 .and. lcplxp .eq. 0)
     .call rx('mkpot: incompatible lso,lcplxp')
      if (isum(nbas,lldau,1) .ne. 0 .and. lcplxp .eq. 0)
     .call rx('mkpot: incompatible ldau,lcplxp')
      stdo = lgunit(1)
      stdl = lgunit(2)
      nsp  = nglob('nsp')
      lxcfun = nglob('lxcf')
      pi   = 4d0*datan(1d0)

      i_copy_size=size(slat%nabc)
      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      vol=slat%vol


      i_copy_size=size(sham%eterms)
      call dcopy(i_copy_size,sham%eterms,1,eterms,1)


C --- Ensure density is positive definite ---
      if (cmdopt('--rhopos',8,0,outs)) then
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              if (dble(smrho(i1,i2,i3,1)) .lt. 0) smrho(i1,i2,i3,1) = 0
              if (dble(smrho(i1,i2,i3,nsp)) .lt. 0) smrho(i1,i2,i3,nsp) = 0
            enddo
          enddo
        enddo
      endif

C --- Printout for smooth background charge ---
      if (qbg .ne. 0) then
        rhobg = (3d0/4d0/pi*vol)**(1d0/3d0)
        call info5(20,1,0,' Energy for background charge'//
     .  ' q=%d, radius r=%;3d :  E = 9/5*q*q/r = %;4d',
     .  qbg,rhobg,1.8d0*qbg*qbg/rhobg,0,0)
      endif

C --- Smooth electrostatic potential ---
Cchp3       call rhomom(nbas,ssite,sspec,orhoat,qmom,vsum)
       call rhomom ( nbas , ssite , sspec , sv_p_orhoat , qmom , vsum 
     . ) 

      vconst = -vsum/vol
C      if (ipr .ge. 40) write (stdo,334) vsum,vconst
C  334 format(' vsum=',f14.6,'   vconst=',f14.6)
      allocate(hpot0_rv(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) hpot0_rv(:)=0.0d0


      i = 1
      if (cmdopt('--oldvc',7,0,outs)) i = 0

cccccccccccccccccccccc
c      smpot(:,:,:,1)=0d0
c      allocate(smpotbk(k1,k2,k3))
c      smpotbk=smpot(:,:,:,1)
ccccccccccccccccccccc

      call smves ( i , nbas , ssite , sspec , slat , k1 , k2 , k3 , 
     .qmom , gpot0 , vval , hpot0_rv , sgp0 , smrho , smpot , vconst 
     ., smq , qsmc , fes , rhvsm , zvnsm , zsum , vesrmt , qbg )


cccccccccccccccccc
c      print *,'vvv fix smves'
c      open(2021,file='smpot.bk')
c      do i3=1,k3
c      do i2=1,k2
c      do i1=1,k1
c        read(2021,"(2d23.15)")smpot(i1,i2,i3,1)
cc        write(2021,"(2d23.15)")smpot(i1,i2,i3,1)
c      enddo
c      enddo
c      enddo
c      close(2021)
c      stop 'eeeeeeeeeeeeeeeeeee'
cccccccccccccccccc

c      if(itest==1) then
c      print *,'vvv: skip smves for smpot'
c      smpot(:,:,:,1)=smpotbk
c      deallocate(smpotbk)
c      stop 'xxxxxxxxxxxxxx'
c      endif
ccccccccccccccccc

      smag = 0
      if (nsp .eq. 2) then
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smag,sum2)
        smag = 2*smag - smq
      endif
      if (allocated(hpot0_rv)) deallocate(hpot0_rv)

C     call zprm('smves',2,smpot,k1,k1,k2)

C --- Add smooth exchange-correlation potential ---
      call dpzero(rveps,2)
      call dpzero(rvvxc,2)
      call dpzero(rvepsv,2)
      call dpzero(rvexv,2)
      call dpzero(rvecv,2)
      call dpzero(rvvxcv,2)
      if (mod(job,1000) .lt. 100) then
        allocate(smvxc_zv(k1*k2*k3*nsp),smvx_zv(k1*k2*k3*nsp),
     &  smvc_zv(k1*k2*k3*nsp),smexc_zv(k1*k2*k3), fxc_rv(3*nbas))
        smvxc_zv(:)=0d0
        smvx_zv(:) =0d0
        smvc_zv(:) =0d0
        smexc_zv(:)=0d0
        fxc_rv(:)  =0d0
C   ... Compute int smrho*exc(smrho) and int smrho*vxc(smrho)
c        call pshpr(ipr-30)
        call smvxc2 ( 0 ,slat, nsp , lxcfun , vol , n1 , n2 , n3 , k1 , k2 
     .  , k3 , smrho , smvxc_zv , smvx_zv , smvc_zv , smexc_zv 
     .  , wdummy , rvepsv , rvexv , rvecv , rvvxcv , vxcavg )
C        call mshdot(vol,nsp,n1,n2,n3,k1,k2,k3,smrho,w(osmvxc),
C     .    tmp,tmp(2))
c        call poppr
        rvepsv(1) = rvepsv(1) + rvepsv(2)
        rvvxcv(1) = rvvxcv(1) + rvvxcv(2)
C   ... Smooth exchange-correlation potential
c        call dpzero ( smvxc_zv , k1 * k2 * k3 * nsp )
c        call dpzero ( smexc_zv , k1 * k2 * k3 )
        smvxc_zv=0d0 
        smexc_zv=0d0 


cccccccccccccccccccccccc
c      open(2021,file='smpot2.bk')
c      do i3=1,k3
c      do i2=1,k2
c      do i1=1,k1
c        read(2021,"(2d23.15)")smpot(i1,i2,i3,1)
c      enddo
c      enddo
c      enddo
c      close(2021)
c      stop 'ccccccccccccccccccc'
ccccccccccccccccccccccc

c      smpot(:,:,:,1)=0d0
c      allocate(smpotbkx(k1,k2,k3))
c      smpotbkx=smpot(:,:,:,1)
cccccccccccccccccc

        call smvxcm ( ssite , sspec , slat , nbas , lfrce , k1 , k2 , 
     .  k3 , smrho , smpot , smvxc_zv , smvx_zv , smvc_zv , 
     .  smexc_zv , repsm , repsmx , repsmc , rmusm , rvmusm , rvepsm 
     .  , fcexc0 , fcex0 , fcec0 , fcvxc0 , fxc_rv )



cccccccccccccccccccccccc
ccc      print *,'vvv readin smpot vxc'
c      open(2021,file='smpot2.bk')
c      allocate(smpotbk(k1,k2,k3))
c      do i3=1,k3
c      do i2=1,k2
c      do i1=1,k1
cc         write(2021,"(2d23.15)")smpot(i1,i2,i3,1)-smpotbkx(i1,i2,i3)
c        read(2021,"(2d23.15)") smpotbk(i1,i2,i3)
c      enddo
c      enddo
c      enddo
c      close(2021)
c      smpot(:,:,:,1) = smpotbkx + smpotbk
c      stop 'xxxxxxx smpot2.bk write'
ccccccccccccccccccccccccc


ccccccccccccc
c        print *,'vvv sum smrho=',sum(smrho(:,:,:,1)),sum(abs(smrho(:,:,:,1)))
c        print *,'vvv repsm=',repsm
cccccccccccc
        if ( lfrce .ne. 0 ) call dpadd ( fes , fxc_rv , 1 , 3 * nbas 
     .  , 1d0 )

        if (allocated(fxc_rv)) deallocate(fxc_rv)
        if (allocated(smexc_zv)) deallocate(smexc_zv)
        if (allocated(smvc_zv)) deallocate(smvc_zv)
        if (allocated(smvx_zv)) deallocate(smvx_zv)
        if (allocated(smvxc_zv)) deallocate(smvxc_zv)

      else
        call dpzero(repsm,2)
        call dpzero(repsmx,2)
        call dpzero(repsmc,2)
        call dpzero(rmusm,2)
        call dpzero(rvmusm,2)
        call dpzero(fcexc0,2)
        call dpzero(fcex0,2)
        call dpzero(fcec0,2)
        call dpzero(fcvxc0,2)
      endif

C --- Make parameters for extended local orbitals ---
      call elocp(nbas,nsp,ssite,sspec,job)

C --- Make local potential at atomic sites and augmentation matrices ---
      rhobg=qbg/vol
Cchp3       call locpot(nbas,nsp,lso,lcplxp,ssite,sspec,slat,orhoat,qmom,vval,
Cchp3      .gpot0,osig,otau,oppi,ppn,hab,vab,sab,vvesat,cpnvsa,repat,
Cchp3      .repatx,repatc,rmuat,rvepva,rvexva,rvecva,rvvxva,rvepsa,rvvxca,
Cchp3      .valvfa,xcore,fcexca,fcexa,fceca,fcvxca,sqloc,sqlocc,saloc,qval,
Cchp3      .qsc,job,rhobg,nlibu,lmaxu,vorb,lldau)
Cchp3        call locpot ( nbas , nsp , lso , lcplxp , ssite , sspec , slat 
Cchp3      . , orhoat , qmom , vval , gpot0 , sv_p_osig , otau , oppi , ppn 
Cchp3      . , hab , vab , sab , vvesat , cpnvsa , repat , repatx , repatc 
Cchp3      . , rmuat , rvepva , rvexva , rvecva , rvvxva , rvepsa , rvvxca 
Cchp3      . , valvfa , xcore , fcexca , fcexa , fceca , fcvxca , sqloc , 
Cchp3      . sqlocc , saloc , qval , qsc , job , rhobg , nlibu , lmaxu , vorb 
Cchp3      . , lldau ) 
Cchp3         call locpot ( nbas , nsp , lso , lcplxp , ssite , sspec , slat 
Cchp3      .  , orhoat , qmom , vval , gpot0 , sv_p_osig , sv_p_otau , oppi 
Cchp3      .  , ppn , hab , vab , sab , vvesat , cpnvsa , repat , repatx , 
Cchp3      .  repatc , rmuat , rvepva , rvexva , rvecva , rvvxva , rvepsa , 
Cchp3      .  rvvxca , valvfa , xcore , fcexca , fcexa , fceca , fcvxca , sqloc 
Cchp3      .  , sqlocc , saloc , qval , qsc , job , rhobg , nlibu , lmaxu , 
Cchp3      .  vorb , lldau ) 
Cchp3          call locpot ( nbas , nsp , lso , lcplxp , ssite , sspec , slat 
Cchp3      .   , orhoat , qmom , vval , gpot0 , sv_p_osig , sv_p_otau , sv_p_oppi 
Cchp3      .   , ppn , hab , vab , sab , vvesat , cpnvsa , repat , repatx , 
Cchp3      .   repatc , rmuat , rvepva , rvexva , rvecva , rvvxva , rvepsa , 
Cchp3      .   rvvxca , valvfa , xcore , fcexca , fcexa , fceca , fcvxca , sqloc 
Cchp3      .   , sqlocc , saloc , qval , qsc , job , rhobg , nlibu , lmaxu , 
Cchp3      .   vorb , lldau ) 
          call locpot ( nbas , nsp , lso , lcplxp , ssite , sspec , slat 
     .    , sv_p_orhoat , qmom , vval , gpot0 , sv_p_osig , sv_p_otau , 
     .    sv_p_oppi , ppn , hab , vab , sab , vvesat , cpnvsa , repat , 
     .    repatx , repatc , rmuat , rvepva , rvexva , rvecva , rvvxva , 
     .    rvepsa , rvvxca , valvfa , xcore , fcexca , fcexa , fceca , fcvxca 
     .    , sqloc , sqlocc , saloc , qval , qsc , job , rhobg , nlibu , 
     .    lmaxu , vorb , lldau ) 





C ... Combine spin-up and spin-down integrals
      repsm(1)  = repsm(1) + repsm(2)
      repsmx(1) = repsmx(1)+ repsmx(2)
      repsmc(1) = repsmc(1)+ repsmc(2)
      rmusm(1)  = rmusm(1) + rmusm(2)
      rvmusm(1) = rvmusm(1) + rvmusm(2)
      fcexc0(1) = fcexc0(1) + fcexc0(2)
      fcex0(1)  = fcex0(1) + fcex0(2)
      fcec0(1)  = fcec0(1) + fcec0(2)
      fcvxc0(1) = fcvxc0(1) + fcvxc0(2)
      repat(1)  = repat(1) + repat(2)
      repatx(1) = repatx(1)+ repatx(2)
      repatc(1) = repatc(1)+ repatc(2)
      rmuat(1)  = rmuat(1) + rmuat(2)
      fcexca(1) = fcexca(1) + fcexca(2)
      fcexa(1)  = fcexa(1) + fcexa(2)
      fceca(1)  = fceca(1) + fceca(2)
      fcvxca(1) = fcvxca(1) + fcvxca(2)

C --- Total energy terms ---
C ... Integral of valence density times xc potential(valence density)
      if (mod(job,1000) .lt. 100) then
        rvepsv(1) = rvepsv(1) + rvepva
        rvexv(1)  = rvexv(1) + rvexva
        rvecv(1)  = rvecv(1) + rvecva
        rvvxcv(1) = rvvxcv(1) + rvvxva
      endif
      if (mod(mod(job/1000,10),2) .eq. 1) then
        rveps(1) = rveps(1) + rvepsa
        rvvxc(1) = rvvxc(1) + rvvxca
      endif


C ... Integral of valence density times estatic potential
      valves = rhvsm + vvesat

C ... Valence density times veff.
C    *Associate term (n0~-n0) Ves(n0~) with local part
C     because of the ppi matrix elements
C    *Also add fcvxc0(1) to smooth part because
C     rvmusm+fcvxc0 is perturbative approximation for rvmusm
C     when cores are not treated perturbatively.
      valfsm = rhvsm + rvmusm(1) - sgp0 - vconst*qbg
ccccccccccccccc
c      write(6,"('vvv: mkpot : rhvsm, sgp0, rhvsm-sgp0',10d13.6)") rhvsm, sgp0, rhvsm-sgp0
ccccccccccccccc
      valftr = valvfa + sgp0
      valfsm = valfsm + fcvxc0(1)
      valftr = valftr - fcvxca(1)
      valvef = valfsm + valftr
cccccccccccccccccc
c      print *,'vvv takao valvef valvfsm valvtr valvfa',valvef,valfsm,valftr,valvfa
cccccccccccccccccc

C ... Integral of core+nucleus times estatic potential
      cpnves = zvnsm + cpnvsa

C ... Total xc energy and potential integral
      focexc = fcexc0(1) - fcexca(1)
C     focex  = fcex0(1)  - fcexa(1)
C     focec  = fcec0(1)  - fceca(1)
      focvxc = fcvxc0(1) - fcvxca(1)
      if (ipr.ge.30 .and. dabs(focexc).gt.1d-6)
     .write (stdo,850) focexc,focvxc
  850 format(' foca xc integrals for spillout charge:',2f12.6)

      repsm(1) = repsm(1) + fcexc0(1)
      repsmx(1)= repsmx(1)+ fcex0(1)
      repsmc(1)= repsmc(1)+ fcec0(1)
      repat(1) = repat(1) - fcexca(1)
      repatx(1)= repatx(1)- fcexa(1)
      repatc(1)= repatc(1)- fceca(1)
      rmusm(1) = rmusm(1) + fcvxc0(1) + fcexc0(1)
      rmuat(1) = rmuat(1) - fcvxca(1) - fcexca(1)

      rhoexc = repsm(1) + repat(1)
      rhoex  = repsmx(1)+ repatx(1)
      rhoec  = repsmc(1)+ repatc(1)
      rhovxc = rmusm(1) + rmuat(1)

C ... Total electrostatic energy
      usm = 0.5d0*(rhvsm+zvnsm)
      uat = 0.5d0*(vvesat+cpnvsa)
      utot = usm + uat
      dq = smq+sqloc + qsmc+sqlocc + qbg -zsum
      amom = smag+saloc

C --- Printout ---
      if (ipr .ge. 20) write(stdo,'(1x)')
      if (ipr .ge. 30) then
        write (stdo,681)
        write (stdo,680) 'rhoval*vef ',valfsm,valftr,valvef,
     .  'rhoval*ves ',rhvsm,vvesat,valves,
     .  'psnuc*ves  ',zvnsm,cpnvsa,cpnves,
     .  'utot       ',usm,uat,utot,
     .  'rho*exc    ',repsm(1),repat(1),rhoexc,
     .  'rho*vxc    ',rmusm(1),rmuat(1),rhovxc,
     .  'valence chg',smq,sqloc,smq+sqloc
        if (nsp .eq. 2)
     .  write (stdo,680) 'valence mag',smag,saloc,amom
        write (stdo,680) 'core charge',qsmc,sqlocc,qsmc+sqlocc
        write (stdo,670) smq+sqloc,qsmc+sqlocc,-zsum,qbg,dq
      endif
  680 format(3x,a,4x,3f17.6)
  681 format(' Energy terms:',13x,'smooth',11x,'local',11x,'total')
  670 format(/' Charges:  valence',f12.5,'   cores',f12.5,
     .'   nucleii',f12.5/'    hom background',f12.5,
     .'   deviation from neutrality: ',f12.5)

ccccccccccccccccccc takao sep2009
c      print *,'vvvvvv ipl=',ipl
c      ipl=1
ccccccccccccccccccc
      if (ipl .ge. 1) then
c         print *,'stdl=',stdl
        write (stdl,710) smq+sqloc,smq,sqloc,qbg,dq
  710   format('fp qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .  '  qbg',f11.6,' dQ',f11.6)
        if (nsp .eq. 2) write (stdl,711) smag+saloc,smag,saloc
  711   format('fp mag',f11.5,'  sm',f11.5,'  loc',f11.5)
        write (stdl,720) rhovxc,rhoexc,utot
  720   format('fp pot  rvxc',f18.7,'  rexc',f18.7,'  rves',f16.7)
        write (stdl,721) rhoex,rhoec
  721   format('fp pot  rex ',f18.7,'  rec',f19.7)
        write (stdl,722) rvvxcv(1),rvepsv(1)
  722   format('fp pot  rvvxcv',f16.7,'  rvexcv',f16.7)
        write (stdl,723) rvexv(1),rvecv(1)
  723   format('fp pot  rvexv',f17.7,'  rvecv',f17.7)
c

c------------------
C      ifet = 1457
C      open(ifet,file='RoVxcLDA')
C      write(ifet,"(d23.16,a)")  rhoexc, ' ! \int rho exc (Ry)'
C      write(ifet,"(d23.16,a)")  rhoex,  ' ! \int rho ex  '
C      write(ifet,"(d23.16,a)")  rhoec,  ' ! \int rho ec  '
C      write(ifet,"(d23.16,a)")  rhovxc, ' ! \int rho vxc '
C      close(ifet)
c------------------

      endif
      if (dabs(dq) .gt. 1d-3 .and. iprint() .gt. 0) call awrit1(
     .' (warning) system not neutral, dq=%d',' ',80,lgunit(1),dq)

      i_copy_size=size(sham%eterms)
      call dcopy(i_copy_size,eterms,1,sham%eterms,1)


      call tcx('mkpot')

      end subroutine mkpot


