module m_locpot
  use m_lmfinit,only: lmxax,nsp,nbas
  use m_MPItk,only: master_mpi
  use m_lgunit,only:stdo
  use m_vxcatom,only: vxcnsp
  public locpot
  real(8),allocatable,public :: rotp(:,:,:,:,:) 
  private
contains
  !- Make the potential axt the atomic sites and augmentation matrices.
  subroutine locpot(orhoat,qmom,vval,gpot0,job,rhobg,nlibu,lmaxu,vorb,lldau,novxc, &!,idipole 
       osig , otau, oppi, ohsozz,ohsopm, ppnl, hab, vab, sab, & 
       vvesat,cpnvsa, rhoexc,rhoex,rhoec,rhovxc, valvef, xcore, sqloc,sqlocc,saloc, qval,qsc )
    use m_density,only: v0pot,v1pot   !output
    use m_density,only: pnzall,pnuall !output
    use m_lmfinit,only:nkaph,lxcf,lhh,nkapii,nkaphh
    use m_lmfinit,only:n0,nppn,nab,nrmx,nkap0,nlmx,nbas,nsp,lso,ispec, sspec=>v_sspec,mxcst4,lmxax
    use m_lmfinit,only:slabl,idu,coreh,ham_frzwf,rsma,alat,v0fix
    use m_uspecb,only:uspecb
    use m_ftox
    use m_struc_def
    use m_augmat,only: augmat
    use m_hansr,only:corprm
    implicit none
    intent(in)::    orhoat,qmom,vval,gpot0,job,rhobg,nlibu,lmaxu,vorb,lldau,novxc
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   orhoat:vector of offsets containing site density
    !i   qmom  :multipole moments of on-site densities (rhomom.f)
    !i   vval  :electrostatic potential at MT boundary; needed
    !i         :to computed matrix elements of local orbitals.
    !i   gpot0 :integrals of local gaussians * phi0~
    !i         :phi0~ is the estatic potential of the interstitial
    !i   job   :1s  digit
    !i         : 1 make core and augmentation matrices
    !i  rhobg  :compensating background density
    !i  lso    :if nonzero, calculate LzSz matrix elements
    !i  lcplxp=1 only now :0 if ppi is real; 1 if ppi is complex
    !i  nlibu  : max number of lda+u blocks
    !i  lmaxu  : max l for U
    !i  vorb   : orbital dependent potential matrices
    !i  lldau  :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat beginning at dmats(*,lldau(ib))
    !o Outputs
    !o   ... the following are summed over all spheres
    !o   vvesat:integral of valence density times electrostatic potential
    !o   cpnvsa:integral of core+nucleus times electrostatic potential
    !o          = int [ (rhoc-z) ves1~ - (rhocsm-rhonsm) ves2~ ]
    !o          where ves1~ is the estat potential of the true density
    !o          where ves2~ is the estat potential of the smooth density
    !o   rhoexc:integral of density times xc energy density
    !o   rhoex :integral of density times exch. energy density
    !o   rhoec :integral of density times corr. energy density
    !o   rhovxc:integral of density times xc potential
    !!o   rvepsv:integral of valence density times exc(valence density)
    !!o   rvexv :integral of valence density times ex(valence density)
    !!o   rvecv :integral of valence density times ec(valence density)
    !!o   rvvxcv:integral of valence density times vxc(valence density)
    !!o   rveps :int rhov*exc(rhotot) (only made if 1000's digit job set)
    !!o   rvvxc :int rhov*vxc(rhotot) (only made if 1000's digit job set)
    !o   valvef:integral (rho1*(v1-2Z/r) - rho2*vsm) ??
    !o   xcore :integral rhoc*(v1-2Z/r)
    !o   sqloc :total valence charge rho1-rho2 in sphere
    !o   saloc :total valence magnetic moment in sphere
    !o         :+core moment for sites with core hole
    !o   sqlocc:total core charge qcor1-qcor2 in sphere
    !o   qval  :nominal total valence charge qv-z
    !o   qsc   :semicore charge from local orbitals
    !o   osig  :augmentation overlap integrals; see augmat
    !o   otau  :augmentation kinetic energy integrals; see augmat
    !o   oppi  :augmentation kinetic + potential integrals; see augmat
    !o   ppnl  :NMTO-like potential parameters
    !o   hab   :integrals of the ham. with true w.f.  See Remarks in augmat
    !o   vab   :integrals of the pot. with true w.f.  See Remarks in augmat
    !o   sab   :integrals of    unity with true w.f.  See Remarks in augmat
    !l Local variables
    !l   lfltwf:T  update potential used to define basis
    !l   iblu  :index to current LDA+U block
    !i   idu   :idu(l+1)=1 => this l has a nonlocal U matrix
    !r Remarks
    !u Updates
    !u   10 Sep 08 Added electric field gradient (AxSv)
    !u   02 Jan 06 adds core magnetic moment to saloc
    !u   09 Nov 05 Convert dmat to complex form
    !u   06 Jul 05 Parameters for extended local orbitals are now input
    !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
    !u   12 Jun 05 Optionally write sphere densities to file
    !u   27 Apr 05 Added LDA+U  (Lambrecht)
    !u    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
    !u   15 Jul 04 First implementation of extended local orbitals
    !u   14 Mar 04 Makes rhov*exc and rhov*vxc
    !u   14 Mar 04 Makes rhov*exc and rhov*vxc
    !u   19 Sep 03 (ATP) Enabled partial core occupation
    !u   02 Oct 02 (WRL) Added background potential
    !u    9 May 02 Added species-specific freezing of potential
    !u    8 May 02 Added rhoex and rhoec (T. Miyake)
    !u    7 May 02 Added rvexv and rvecv (T. Miyake)
    !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
    !u   17 Sep 01 Returns qsc.  Altered argument list.
    !u   28 Aug 01 Extended to local orbitals.  Altered argument list.
    !u   15 Aug 01 Generates rvepsv and rvvxcv.  New argument list
    !u   18 Apr 01 Added ability to exclude exchange-correlation potential
    !u   20 Feb 01 Added ppnl to potential parameters generated
    !u   13 Jun 00 spin polarized
    !u    1 May 00 Adapted from nfp locpot.f
    !u   17 Jun 98 Adapted from DLN to run parallel on SGI
    ! ----------------------------------------------------------------------
    integer::  job,ibx,ir,isp,l,lm
    type(s_rv1) :: orhoat(3,nbas)
    type(s_cv1) :: oppi(3,nbas)
    type(s_sblock):: ohsozz(3,nbas),ohsopm(3,nbas)
    type(s_rv1) :: otau(3,nbas)
    type(s_rv1) :: osig(3,nbas)
    integer :: nlibu,lmaxu,lldau(nbas),iblu
    double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    real(8):: qmom(1) , vval(1)
    double precision :: cpnvsa,rhoexc(nsp),rhoex(nsp),rhoec(nsp),rhovxc(nsp), &
         qval,sqloc,sqlocc,saloc, & !,focvxc(nsp)focexc(nsp),focex(nsp),focec(nsp),
         valvef,vvesat,xcore,& !rvvxc, & !,rveps rvepsv, ,rvexv,rvecv,rvvxcv
         hab(3,3,n0,nsp,nbas),vab(3,3,n0,nsp,nbas),sab(3,3,n0,nsp,nbas), &
         gpot0(1),ppnl(nppn,n0,nsp,nbas),rhobg
    character spid*8
    integer :: lh(nkap0),nkapi,nkape,k
    double precision :: eh(n0,nkap0),rsmh(n0,nkap0)
    double precision :: ehl(n0),rsml(n0)
    double precision :: rofi(nrmx),rwgt(nrmx), &
         gpotb(81),z,a,rmt,qc,ceh,rfoc, &
         qcorg,qcorh,qsc,cofg,cofh,qsca,rg,qv,cpnvs, &
         qloc,qlocc,xcor, aloc,alocc!,rvexl, rvecl,rvvxvl,rvvxtl !,rvepvl,rveptl
    real(8),pointer:: pnu(:,:),pnz(:,:)
    ! ... for sm. Hankel tails
    double precision :: rs3,vmtz
    character chole*8
    integer :: kcor,lcor
    double precision :: qcor(2),qc0,qsc0
    real(8),allocatable:: wk(:),rhol1(:),rhol2(:), v1(:),v2(:),v1es(:),v2es(:),efg(:,:),zz(:)
    logical :: lfltwf
    integer:: i,nglob,ipr,iprint,j1,ib,is,lmxl,lmxa,nr,lmxb,kmax,lfoc,nrml,nlml,ifivesint,ifi
    logical,save:: secondcall=.false.
    logical :: phispinsym,cmdopt0,readov0,v0write,novxc
    integer::lsox
    character*20::strib
    character strn*120
    real(8):: ov0mean,pmean
    real(8),parameter:: pi = 4d0*datan(1d0),srfpi = dsqrt(4d0*pi),y0 = 1d0/srfpi  !integer:: idipole
    call tcn('locpot')
    ipr = iprint()
    if (ipr >= 30) write(stdo,"('  locpot:')")
    k = nrmx*nlmx*nsp
    allocate(rhol1(k),rhol2(k),v1(k),v2(k),v1es(k),v2es(k),efg(5,nbas),zz(nbas))
    xcore   = 0d0
    if(master_mpi) open(newunit=ifivesint,file='vesintloc',form='formatted',status='unknown')
    if(allocated(rotp)) deallocate(rotp)
    allocate(rotp(0:lmxax,nsp,2,2,nbas))
    ibblock: block
      real(8):: valvs(nbas),cpnvs(nbas),valvt(nbas)
      real(8)::qloc(nbas),aloc(nbas),qlocc(nbas),alocc(nbas) 
      real(8):: rhexc(nsp,nbas),rhex(nsp,nbas),rhec(nsp,nbas),rhvxc(nsp,nbas)
      real(8):: xcor(nbas),qv(nbas),qsca(nbas)
      valvs=0d0;cpnvs=0d0;valvt=0d0 
      qloc=0d0;aloc=0d0;qlocc=0d0;alocc=0d0
      rhexc=0d0;rhex=0d0;rhec=0d0;rhvxc=0d0 ;xcor=0d0;qv=0d0;qsca=0d0
      iblu = 0
      j1 = 1
      ibloop: do  ib = 1, nbas
         is=ispec(ib) 
         pnu=>pnuall(:,:,ib)
         pnz=>pnzall(:,:,ib)
         z=sspec(is)%z
         qc=sspec(is)%qc
         rg=sspec(is)%rg
         spid=slabl(is) 
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt
         lmxa=sspec(is)%lmxa
         lmxl=sspec(is)%lmxl
         lmxb=sspec(is)%lmxb
         zz(ib)=z
         if (lmxa == -1) cycle ! floating orbital
         kmax=sspec(is)%kmxt
         lfltwf = (.not.mxcst4(is)).and.(.not.ham_frzwf).and.job==1 ! modify b.c. of Rad.wave func.
         call corprm(is, qcorg,qcorh,qsca(ib),cofg,cofh,ceh,lfoc,rfoc,z)
         chole=coreh(is)
         call gtpcor(is, kcor,lcor,qcor) !qcor(1:2) is meaningful only when kcor/=0 
         call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor, qc0,qv(ib),qsc0)
         if (qsc0 /= qsca(ib) .OR. qc /= qc0-qsc0) then
            if(iprint()>0)write(stdo,ftox)' is=',is,'qsc0=',ftof(qsc0),'qsca',ftof(qsca(ib)),'qc',ftof(qc),'qc0',ftof(qc0)
            call rxs('problem in locpot -- possibly low LMXA or orbital mismatch, species ',spid)
         endif
         nlml = (lmxl+1)**2
         nrml = nr*nlml
         if (ipr >= 20) then
            write(stdo,"('   site',i3,'  z=',f5.1,'  rmt=',f8.5,'  nr=',i3,'   a=',f5.3, &
                 '  nlml=',i2,'  rg=',f5.3,'  Vfloat=',l1)") ib,z,rmt,nr,a,nlml,rg,lfltwf
            if (kcor/= 0) then
               if(sum(abs(qcor)) /= 0 ) then
                  if(ipr>=30)write(stdo,ftox)' core hole: kcor=',kcor,'lcor=',lcor,'qcor amom=',ftof(qcor)
               endif
            endif
         endif
         call rxx(nr .gt. nrmx,  'locpot: increase nrmx')
         call rxx(nlml .gt. nlmx,'locpot: increase nlmx')
         call radmsh(rmt,a,nr,rofi)
         call radwgt(rmt,a,nr,rwgt)
         !   ... Write true density to file rhoMT.ib
         if(cmdopt0('--wrhomt'))call wrhomt('rhoMT.','density',ib,orhoat(1,ib)%v,rofi,nr,nlml,nsp)
         call locpt2(z,rmt,rg,a,nr,nsp,cofg,cofh & ! Make potential and energy terms at this site ---
              ,ceh,rfoc,lfoc,nlml,qmom ( j1 ),vval ( j1 ),rofi &
              ,rwgt,orhoat( 1,ib )%v,orhoat( 2,ib )%v,&
              orhoat( 3,ib )%v,rhol1,rhol2,v1,v2,v1es,v2es, &
              valvs(ib),cpnvs(ib),rhexc(:,ib),rhex(:,ib),rhec(:,ib),rhvxc(:,ib),&
              valvt(ib),xcor(ib) ,qloc(ib),qlocc(ib),aloc(ib),alocc(ib),gpotb,& 
              rhobg,efg(1,ib),ifivesint,lxcf) 
         !! write density 1st(true) component and counter components.
         if(cmdopt0('--density') .AND. master_mpi .AND. secondcall) then
            write(stdo,"(' TotalValenceChange diff in MT; ib,\int(rho2-rho1)=',i5,f13.5)") ib,qloc(ib)
            write(strib,'(i10)') ib
            open(newunit=ibx,file='rho1MT.'//trim(adjustl(strib)))
            do isp=1,nsp
               do ir=1,nr
                  write(ibx,"(i2,d13.6,256f13.6)")isp,rofi(ir),(rhol1(ir+(lm-1)*nr+(isp-1)*nr*nlml),lm=1,nlml)
               enddo
            enddo
            close(ibx)
            open(newunit=ibx,file='rho2MT.'//trim(adjustl(strib)))
            do isp=1,nsp
               do ir=1,nr
                  write(ibx,"(i2,d13.6,256f13.6)")isp,rofi(ir),(rhol2(ir+(lm-1)*nr+(isp-1)*nr*nlml),lm=1,nlml)
               enddo
            enddo
            close(ibx)
         endif
         
         !   ... Write true potential to file vtrue.ib
         if(cmdopt0('--wpotmt'))call wrhomt('vtrue.','potential',ib,v1,rofi,nr,nlml,nsp)
         if (lfltwf) then !   ... Update the potential used to define basis set
            do i = 0, nsp-1
               v0pot(ib)%v(1+nr*i: nr+nr*i) = y0*v1(1+nr*nlml*i : nr+nr*nlml*i)
            enddo
         endif
         
         phispinsymB: block ! spin averaged oV0 to generate phi and phidot. takaoAug2019
           phispinsym= cmdopt0('--phispinsym')
           if(phispinsym) then
              if(master_mpi .AND. nsp==2) then
                write(6,*) 'locpot: --phispinsym mode: use spin-averaged potential for phi and phidot'
              endif
              do ir=1,nr
                 ov0mean = 0d0
                 do isp=1,nsp
                    ov0mean = ov0mean + v0pot(ib)%v( ir + nr*(isp-1) )
                 enddo
                 ov0mean = ov0mean/nsp
                 do isp=1,nsp
                    v0pot(ib)%v(ir + nr*(isp-1))= ov0mean
                 enddo
              enddo
           endif
         endblock phispinsymB

         v0fixblock: block ! experimental case --v0fix
           character charext*8
           ! v0write=cmdopt0('--v0write')
           ! if(v0write) then !v0write is rdovfa 
           !    do ir=1,nr
           !       ov0mean = 0d0
           !       do isp=1,nsp
           !          ov0mean = ov0mean + v0pot(ib)%v( ir + nr*(isp-1) )
           !       enddo
           !       ov0mean = ov0mean/nsp !spin averaged
           !       do isp=1,nsp
           !          v0pot(ib)%v(ir + nr*(isp-1))= ov0mean
           !       enddo
           !    enddo
           !    open(newunit=ifi,file='v0pot.'//trim(charext(ib)),form='unformatted')
           !    write(ifi) v0pot(ib)%v(1:nr)
           !    close(ifi)
           !    if(ib==nbas) readov0= .TRUE.
           !    call rx0('end of v0write')
           ! endif
           if(v0fix) then
              inquire(file='v0pot.'//trim(charext(ib)),exist=readov0)
              write(6,*)'v0fixmode=',readov0,ib,nr
              if(readov0) then
                 v0potb:block
                   real(8):: ov0(nr)
                   open(newunit=ifi,file='v0pot.'//trim(charext(ib)),form='unformatted')
                   read(ifi) ov0(1:nr)
                   close(ifi)
                   do ir=1,nr
                      do isp=1,nsp
                         v0pot(ib)%v(ir + nr*(isp-1))= ov0(ir)
                      enddo
                   enddo
                 endblock v0potb
              else
                 call rx('no v0pot files')
              endif
           endif
         endblock v0fixblock

         if(master_mpi .AND. nsp==2)then
            do l=0,lmxa
               write(6,"(' ibas l=',2i3,' pnu(1:nsp) pnz(1:nsp)=',4f10.5)") ib,l,pnu(l+1,1:nsp),pnz(l+1,1:nsp)
            enddo
         endif
         
         do  i = 0, nsp-1 ! Store the potential used in mkrout to calculate the core
            v1pot(ib)%v(1+nr*i: nr+nr*i) = y0*v1(1+nr*nlml*i : nr+nr*nlml*i)
         enddo
         if (lfoc==0) xcore = xcore + xcor(ib)
         if (kcor/=0) then !  Check for core moment mismatch ; add to total moment
            if (dabs(qcor(2)-alocc(ib)) > 0.01d0) then
               if(ipr>=10) write(stdo,ftox) ' (warning) core moment mismatch spec=',is,&
                    'input file=',qcor(2),'core density=',alocc
            endif
         endif
         !   --- Make augmentation matrices sig, tau, ppi ---
         if (job==1) then !     ... Smooth Hankel tails for local orbitals
            rsmh= 0d0
            eh  = 0d0
            call uspecb(is,rsmh,eh)
            nkapi=nkapii(is)
            nkape=nkaphh(is)
            if (nkape > nkapi) then
               rsml=rsmh(:,nkaph)
               ehl =  eh(:,nkaph)
            endif
            if (novxc) then !  ... Use effective potentials with modified xc
               v1=v1es 
               v2=v2es 
            endif
            lsox=lso
            !          if(idipole/=0) then
            !             call adddipole(v1,rofi,nr,nlml,nsp,idipole,ssite(ib)%pos(idipole)*alat)
            !             call adddipole(v2,rofi,nr,nlml,nsp,idipole,ssite(ib)%pos(idipole)*alat)
            !          endif
            if( .NOT. novxc .AND. cmdopt0('--socmatrix') ) lsox=1
            if (ipr >= 20) write(stdo,467) y0*(gpot0(j1)-gpotb(1))
467         format('     potential shift to crystal energy zero:',f12.6)
            call augmat ( z,rmt,rsma(is),lmxa,pnu,pnz,kmax,nlml, a,nr,nsp,lsox,rwgt,& !rofi,
                 v0pot(ib)%v,v1,v2,gpotb,gpot0 ( j1 ),nkaph,nkapi,&
                 lmxb,lhh(:,is),eh,rsmh, ehl,rsml,rs3,vmtz, lmaxu, vorb, lldau(ib), idu, &
                 iblu, &
                 osig(1,ib), otau(1,ib), oppi(1,ib), ohsozz(1,ib), ohsopm(1,ib), &
                 ppnl(1,1,1,ib), hab(1,1,1,1,ib),vab (1,1,1,1,ib), sab(1,1,1,1,ib),rotp(0,1,1,1,ib) )
         endif
         j1 = j1+nlml
      enddo ibloop
      vvesat = sum(valvs)
      cpnvsa = sum(cpnvs)
      valvef=  sum(valvt)
      sqloc  =  sum(qloc)
      saloc  =  sum(aloc)
      if(kcor/=0) saloc = sum(aloc) + qcor(2)
      sqlocc =  sum(qlocc)
      rhoexc =  sum(rhexc,dim=2)
      rhoex  =  sum(rhex,dim=2)
      rhoec  =  sum(rhec,dim=2)
      rhovxc =  sum(rhvxc,dim=2)
      qval = sum(qv)+sum(qsca)
      qsc  = sum(qsca)
    endblock ibblock
    if(cmdopt0('--density') .AND. master_mpi) secondcall= .TRUE. 
    if(master_mpi) close(ifivesint)
    if (ipr > 40) call elfigr(nbas,stdo,zz,efg)!  Electric field gradient
    deallocate(efg,zz,rhol1,rhol2,v1,v2,v1es,v2es)
    call tcx('locpot')
  end subroutine locpot

  subroutine locpt2(z,rmt,rg,a,nr,nsp,cofg,cofh,ceh,rfoc,lfoc, &
       nlml,qmom,vval,rofi,rwgt,rho1,rho2,rhoc,rhol1,rhol2,v1,v2,v1es, &
       v2es,vvesat,cpnves,rhoexc,rhoex,rhoec,rhovxc, valvef,xcore,qloc, & 
       qlocc,aloc,alocc,gpotb,rhobg,efg,ifivesint,lxcfun) 
    use m_hansr,only:hansmr
    use m_ftox
    !- Makes the potential at one site, and associated energy terms.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   z     :nuclear charge
    !i   rmt   :augmentation radius
    !i   rg    :smoothing radius for compensating gaussians used to
    !i         :correct the multipole moments of local smooth density
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
    !i         :cofg = Y0 * qcorg
    !i   cofh  :coefficient to Hankel part of pseudocore density (corprm)
    !i   ceh   :energy of hankel function to fit core tail
    !i   rfoc  :smoothing radius for hankel head fitted to core tail
    !i   lfoc  :switch specifying treatment of core density.
    !i         :0 => core w.f. have val,slo = 0 at rmt
    !i         :1 => core included explicitly with valence
    !i         :2 => core included perturbatively
    !i   nlml  :L-cutoff for rho1,rho2
    !i   qmom  :multipole moments of on-site densities (rhomom.f)
    !i   vval  :boundary condition for estat potential at MT boundary.
    !i         :See remarks.
    !i   rofi  :radial mesh points for density and potential
    !i   rwgt  :radial mesh weights to integrate function on radial mesh
    !i         :Integral f(r) dr = sum_i f_i * wt_i
    !i   rho1  :local valence density = sum_ilm rho1(ilm) * r**2  Y_L(ilm),
    !i         :on radial mesh.
    !i   rho2  :local smoothed valence density, defined as rho1
    !i         :Local atomic valence density is rho1-rho2
    !i   rhoc  :core density times 4*pi*r*r
    !i   rhobg: compensating background density
    !o Outputs
    !o   cpnves:integral of core+nucleus times electrostatic potential
    !o          = int [ (rhoc-z) ves1~ - (rhocsm-rhonsm) ves2~ ]
    !o          where ves1~ is the estat potential of the true density
    !o          where ves2~ is the estat potential of the smooth density
    !o   vvesat:integral of valence density times electrostatic potential
    !o   rhoexc:integral of valence density times xc energy density
    !o   rhoex :integral of valence density times exch. energy density
    !o   rhoec :integral of valence density times corr. energy density
    !o   rhovxc:integral of valence density times xc potential
    !o   valvef:integral valence density times full potential:
    !o         :used to compute double-counting for kinetic energy.
    !o         : = vefv1 - vefv2; defined in Local variables below
    !o   xcore :integral rhoc*(v1-2Z/r)
    !o   qloc  :total valence charge in sphere rho1 - rho2
    !o   aloc  :total valence magnetic moment in sphere
    !o   qlocc :total core charge in sphere qcor1 - qcor2
    !o   alocc :total core magnetic moment in sphere
    !o   gpotb :integrals of gaussians(radius rg) times smooth ves.
    !o         :Here smooth ves is the electrostatic potential of the
    !o         :compensated sm density rhol2, w/ b.c. ves(rmax)=0
    !o         :This is a local analog of gpot0 generated for the smooth
    !o         :density.  Here gpotb is generated in a Y_lm expansion
    !o         :of the gaussian orbitals
    !o   rhol1 :full true electron density, rho1 + rhoc
    !o   rhol2 :full smooth density, i.e. uncompensated rho2 plus
    !o         :compensating gaussians + pseudocore charge
    !o   v1    :true potential at this site, excluding nuclear contribution
    !o         :but including background contribution; see Remarks
    !o         :See Remarks concerning boundary condition at rmt
    !o   v2    :Total smooth potential including background potential
    !o         : = ves[rhol2] + ...
    !o         :   ... lfoc=0 : vxc(rho2)
    !o         :       lfoc=1 : vxc(rho2+sm-rhoc)
    !o         :Apart from differences in l-truncation in the XC potential
    !o         :(i.e. true and local sm. densities are different),
    !o         :v1(nr,1,1)-2*z/y0/rmt = v2(nr,1,1)
    !o         :See Remarks concerning boundary condition at rmt
    !o  efg    :l=2 potential at nucleus
    !l Local variables
    !l   rhol1 :full electron density, rho1 + rhoc
    !l   rhol2 :full smooth compensated density rho2 + gval + gcor + gnuc
    !l   rhochs:Part of smooth core density contains Hankel tails
    !l         :Used to make vxc[n2] and exc[n2]
    !l         :rhochs = srfpi*cofh*xi(0)*r*r
    !l   vefv1 :int (valence true density)  * (total true potential):
    !l         : = int n1 * v1~
    !l   vefv2 :int (valence smooth density)* (total smooth potential):
    !l         : = int (n2*v2~) +  sum_L qmom_L gpotb_L
    !l   qcor1 :true core charge inside rmt
    !l   qcor2 :smoothed core charge inside rmt
    !l   rhonsm:nuclear density smoothed into gaussian of width rg
    !l   rhocsm:core density smoothed into gaussian of width rg
    !r Remarks
    !r   On boundary conditions for the estat potential vval at the MT
    !r   boundary.  In the original formulation (see Springer book) the
    !r   the total hamiltonian and energy did not depend on the choice of
    !r   vval, as the contribution from true and smooth potentials exactly
    !r   cancel.  However, this is no longer true when local orbitals are
    !r   present that have no corresponding smooth part.
    !r
    !r   The potential due to uniform background density rhobg is added
    !r                 vbg= -(4pi/3)*(rmt^2-r^2)*rhobg to v1 and v2
    !r   in spherical components .
    !u Updates
    !u   30 Sep 08 Return l=2 potential at nucleus in efg (A Svane)
    !u   12 Aug 04 First implementation of extended local orbitals
    !u   19 Sep 02 added uniform bkg potential interactions (WRL)
    !u    8 Feb 02 rhoex and rhoec (T. Miyake)
    ! ----------------------------------------------------------------------
    implicit none
    integer :: nr,nsp,lfoc,nlml
    double precision :: z,rmt,rg,a,cofg,cofh,ceh,rfoc,xcore,qloc,qlocc, &
         aloc,alocc,rhoexc(nsp),rhoex(nsp),rhoec(nsp),rhovxc(nsp), &
         valvef,vvesat, cpnves,rhobg
    real(8):: rofi(nr),rwgt(nr), &
         qmom(nlml),vval(nlml),gpotb(nlml), &
         rho1(nr,nlml,nsp),rhol1(nr,nlml,nsp), &
         rho2(nr,nlml,nsp),rhol2(nr,nlml,nsp), &
         v1(nr,nlml,nsp),v1es(nr,nlml,nsp), &
         v2(nr,nlml,nsp),v2es(nr,nlml,nsp), &
         wk(nr,nlml,nsp),rhoc(nr,nsp)
    double precision :: efg(5)
    integer :: ipr,iprint,ll,i,isp,ilm,l,lxcfun,nglob,nrml
    double precision :: rhochs(nr*2),rhonsm(nr),df(0:20),cof(nlml), &
         rhocsm(nr),tmp(nsp),xil(0:0),xill(nr) !xi(0:20,2),
    double precision :: afoc,ag,b,cof0,fac,qv1,qv2,qcor1,qcor2, &
         r,rep1(nsp),rep2(nsp),rep1x(nsp),rep2x(nsp),rep1c(nsp),rep2c(nsp), &
         rhves1,rhves2,rmu1(nsp),rmu2(nsp),rvs1,rvs2, &
         rvsm(nsp),rvtr(nsp),samh,sfac,sgpotb,sum1,sum2,sumg,sumh,top, &
         ves1,vales1,vales2,vcpn1,vcpn2,vefc1,vefv1,vefv2,vesc1,vesc2, &
         vesn1,vesn2,vnucl,vsum,vtr,ddot,a1,a2,smrhoc
    double precision :: qs(nsp)
    real(8),parameter:: pi = 4d0*datan(1d0),srfpi = dsqrt(4d0*pi),y0 = 1d0/srfpi
    real(8):: ves1int,ves2int, w2(nsp),fl(1,1,1),gnu
    logical:: debug=.false.,topl
    integer:: ifivesint
    real(8):: gg(nr)
    call tcn('locpt2')
    call stdfac(20,df)
    ipr = iprint()
    alocc = 0d0
    b = rmt/(dexp(a*nr-a)-1)
    ag = 1d0/rg
    afoc = 1d0/rfoc
    nrml = nr*nlml
    fac = 4d0*pi*(ag*ag/pi)**1.5d0
    gg =  rofi**2 * dexp(-ag*ag*rofi**2) ! ... Renormalize gaussian
    sumg=fac*sum(rwgt(2:nr)*gg(2:nr))
    if(dabs(sumg-1d0)>1d-4)write(stdo,ftox)' locpot (warning): large gaussian, integral=',ftod(sumg)
    sfac = 1d0/sumg
    fac = fac*sfac
    sumh = 0d0
    !  Smooth nuc. and core rho, sm Hankel portion, true & smooth core q
    do i=1,nr
       call hansmr(rofi(i),ceh,afoc,xil,0)
       xill(i)=xil(0)
    enddo
    ! --- Make core and nucleus pseudodensities ---
    rhonsm(:) = -z*fac *gg(:) ! nucleus Gaussian (negative sign)
    rhochs(:) = srfpi*cofh*xill(:)*rofi(:)**2 !pseudocore
    rhocsm(:) = srfpi*cofg*fac *gg(:)   + rhochs(:) ! pcore= pseudocore - Gaussian
    sumh  = sum(rwgt*rhochs)
    samh = -y0*cofh*4d0*pi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
    if(ipr>=20.AND.dabs(samh)>1d-6) write(stdo,ftox)'    sm core charge in MT=',ftof(sumh),&
         '=total-spillout=',ftof(samh),'-',ftof(samh-sumh)
    rhol1 = rho1
    rhol1(:,1,:)= rhol1(:,1,:)+ y0*rhoc(:,:) ! full electron density = rho1 + rhoc 
    ! --- rhol2 = full smooth compensated density rho2+gval +gnuc + pcore 
    !     gval : qmom * compensating gaussians
    !     Distribute core+nuclear charge equally over spins
    do  isp = 1, nsp
       do  ilm = 1, nlml
          l = ll(ilm)
          cof(ilm) = qmom(ilm)*4d0*pi/df(2*l+1)
          fac = sfac*(ag*ag/pi)**1.5d0 * (2d0*ag*ag)**l
          rhol2(:,ilm,isp) = rho2(:,ilm,isp) + cof(ilm)*fac* rofi(:)**l*gg(:)/nsp
       enddo
       rhol2(:,1,isp) = rhol2(:,1,isp) + y0/nsp*(rhocsm(:)+rhonsm(:))
    enddo
    if(nsp==2) then !magnetic moments 
       a1    = srfpi*ddot(nr,rwgt,1,rho1(:,1,1)-rho1(:,1,2),1)
       a2    = srfpi*ddot(nr,rwgt,1,rho2(:,1,1)-rho2(:,1,2),1)
       aloc  = a1-a2
       alocc = ddot(nr,rwgt,1,rhoc(:,1)-rhoc(:,2),1)
    else
       aloc=0d0
       alocc=0d0
    endif
    rhototal: block
      real(8):: ddot,rho1t(nr,nlml),rhol1t(nr,nlml), rho2t(nr,nlml),rhol2t(nr,nlml),rhoct(nr)
      rho1t =sum(rho1(:,:,1:nsp),dim=3) !spin sum total density
      rho2t =sum(rho2(:,:,1:nsp),dim=3)
      rhol1t=sum(rhol1(:,:,1:nsp),dim=3)
      rhol2t=sum(rhol2(:,:,1:nsp),dim=3)
      rhoct =sum(rhoc(:,1:nsp),dim=2)
      rhol1t(:,1)=rhol1t(:,1)+srfpi*rhobg*rofi(:)**2 ! ... Add background density to spherical rhol1 and rhol2
      rhol2t(:,1)=rhol2t(:,1)+srfpi*rhobg*rofi(:)**2
      qv1   = srfpi*ddot(nr,rwgt,1,rho1t,1) ! ... Sphere charges; also check sphere neutrality for safety
      qv2   = srfpi*ddot(nr,rwgt,1,rho2t,1)
      qcor1 =       ddot(nr,rwgt,1,rhoct,1)
      qcor2 =       ddot(nr,rwgt,1,rhocsm,1)
      qlocc = qcor1-qcor2
      qloc  = qv1-qv2
      sum1  = srfpi*ddot(nr,rwgt,1,rhol1t,1) - z
      sum2  = srfpi*ddot(nr,rwgt,1,rhol2t,1)
      qlocc = qcor1-qcor2
      qloc  = qv1-qv2
      if(dabs(sum1-sum2)>1d-6)call rx1('locpt2: sphere not neutral: charge = %d',sum1-sum2)
      !     v1=Ves[rho1t]: true ES pot without nuclear contribution
      call poinsp(z,vval,nlml,a,b,v1,rofi,rhol1t,wk,nr,rvs1,rhves1,  vnucl,vsum)
      if (nlml >= 9 .AND. z > 0.01) then
         efg(1:5)=v1(5,5:9,1)/rofi(5)**2
      else
         efg(1:5)=0d0
      endif
      call poinsp(0d0,vval,nlml,a,b,v2,rofi,rhol2t,wk,nr,rvs2,rhves2, vnucl,vsum)! v2=Ves[rhol2=rho2+gval+gcor+gnuc]
      ! --- gpotb = integrals of compensating gaussians times smooth ves ---
      sgpotb = 0d0
      do  ilm = 1, nlml
         l = ll(ilm)
         cof0 = 4d0*pi/df(2*l+1)
         fac = sfac*(ag*ag/pi)**1.5d0 * (2d0*ag*ag)**l
         gpotb(ilm) = sum(rwgt*v2(:,ilm,1)*cof0*fac*rofi(:)**l* gg(:) )
         sgpotb = sgpotb + qmom(ilm)*gpotb(ilm)
      enddo
      ! --- Electrostatic integrals involving spherical terms only ---
      vesc1 = sum(rwgt(2:nr)*rhoct(2:nr)*(y0*v1(2:nr,1,1) - 2d0*z/rofi(2:nr)))
      vesn2 = sum(rwgt*rhonsm*y0*v2(:,1,1))
      vesc2 = sum(rwgt*rhocsm*y0*v2(:,1,1))
      vnucl = sum(rwgt(2:nr)*rhol1t(2:nr,1)*(1d0/rofi(2:nr)-1d0/rmt))
      ves1int = 4d0*pi*(sum(rwgt*y0*v1(:,1,1)*rofi(:)**2) - z*rofi(nr)**2)
      ves2int = 4d0*pi*sum(rwgt*y0*v2(:,1,1)*rofi(:)**2)
      if(master_mpi)write(ifivesint,"(3f23.15,a)")ves1int-ves2int,ves1int,ves2int,' ! vesint1-vesint2 ves1int ves2int'
      vnucl = 2d0*srfpi*vnucl + 2d0*z/rmt + y0*vval(1)
      vesn1 = -z*vnucl
      vales1 = rvs1-vesc1 ! ... Valence density times electrostatic potential
      vales2 = rvs2-vesn2-vesc2
      vvesat = vales1-vales2
      vcpn1  = vesc1 + vesn1 ! ... Core plus nucleus times estatic potential
      vcpn2  = vesn2 + vesc2
      cpnves = vcpn1 - vcpn2
    endblock rhototal
    if (nsp == 2) then
       v1(:,:,2)=v1(:,:,1) 
       v2(:,:,2)=v2(:,:,1) 
    endif
    v1es=v1 !! ... Preserve potentials without exchange-correlation
    v2es=v2
    call pshpr(max(ipr-31,min(ipr,10)))!  if(debug) write(6,'(a)')' === rho1 valence true density ==='
    call vxcnsp(0,a,rofi,nr,rwgt,nlml,nsp,rho1,lxcfun,rep1,rep1x,rep1c,rmu1,v1,fl,qs)! if(debug) write(6,'(a)')' === rho2 valence counter density ==='
    call vxcnsp(0,a,rofi,nr,rwgt,nlml,nsp,rho2,lxcfun,rep2,rep2x,rep2c,rmu2,v1,fl,qs)
    v1=v1es 
    call poppr
    ! --- Add xc potentials to v1 and v2 ---
    call pshpr(max(ipr-11,min(ipr,10))) !if(debug) write(stdo,'(a)')' === rhol1 valence+core density. rho2->valence+smooth ==='
    call vxcnsp(0,a,rofi,nr,rwgt,nlml,nsp,rhol1,lxcfun,rep1,rep1x,rep1c,rmu1,v1,fl,qs)
    if(lfoc == 0) then
       call vxcnsp(0,a,rofi,nr,rwgt,nlml,nsp,rho2,lxcfun,rep2,rep2x,rep2c,rmu2,v2,fl,qs)
    else if (lfoc == 1) then !frozen core method.
       if (ipr > 40) print *, 'exchange for smooth density, foca=1:'
       do  isp = 1, nsp
          rho2(1:nr,1,isp)=rho2(1:nr,1,isp) + y0/nsp*rhochs(1:nr)
       enddo
       call vxcnsp(0,a,rofi,nr,rwgt,nlml,nsp,rho2,lxcfun,rep2,rep2x,rep2c,rmu2,v2,fl,qs)
       do  isp = 1, nsp
          rho2(1:nr,1,isp)=rho2(1:nr,1,isp)-y0/nsp*rhochs(1:nr)
       enddo
    else
       call rxi('locpt2: cannot handle lfoc = ',lfoc)
    endif
    call poppr 
    vefc1 = sum([(sum(rwgt(2:nr)*rhoc(2:nr,isp)*(y0*v1(2:nr,1,isp)-2d0*z/rofi(2:nr))),isp=1,nsp)])
    vefv2 = 0d0
    vefv1 = 0d0
    if (ipr>=40 .AND. nsp == 1) write(stdo,"(/' ilm',09x,'rho*vtrue',07x,'rho*vsm')")
    if (ipr>=40 .AND. nsp == 2) write(stdo,"(/' ilm',19x,'rho*vtrue',30x,'rho*vsm'/13x, &
         'spin1',7x,'spin2',7x,'tot',11x,'spin1',7x,'spin2',7x,'tot')")
    do  ilm = 1, nlml ! --- Integrals involving the full nonspherical potential ---
       do  isp = 1, nsp
          rvtr(isp)= sum(rwgt(2:nr)*rho1(2:nr,ilm,isp)*v1(2:nr,ilm,isp))
          if(ilm==1) rvtr(isp)=rvtr(isp)- srfpi*2d0*z*sum(rwgt(2:nr)*rho1(2:nr,ilm,isp)/rofi(2:nr))
          rvsm(isp)= sum(rwgt(2:nr)*rho2(2:nr,ilm,isp)*v2(2:nr,ilm,isp))
          vefv1 = vefv1 + rvtr(isp)
          vefv2 = vefv2 + rvsm(isp)
       enddo
       topl = dmax1(dabs(rvsm(1)),dabs(rvtr(1)))>1d-6.and.ipr>=40
       if(topl.AND.nsp == 1)write(stdo,"(i4,3x,2f15.6)") ilm,rvtr(1),rvsm(1)
       if(topl.AND.nsp == 2)write(stdo,"(i4,3x,3f12.6,2x,3f12.6,2x)") ilm,rvtr(1),rvtr(2),&
            rvtr(1)+rvtr(2),rvsm(1),rvsm(2),rvsm(1)+rvsm(2)
    enddo
    rhoexc = rep1 - rep2
    rhoex  = rep1x - rep2x
    rhoec  = rep1c - rep2c
    rhovxc = rmu1 - rmu2
    xcore  = vefc1
    vefv2  = vefv2 + sgpotb
    valvef = vefv1 - vefv2
    if (ipr >= 40) then
       write(stdo,"(/' local terms:     true',11x,'smooth',9x,'local')")
       write(stdo,"(' rhoeps:  ',3f15.6/' rhomu:   ',3f15.6)") sum(rep1),sum(rep2),sum(rhoexc),rmu1(1),rmu2(1),rhovxc(1)
       if(nsp==2)write(stdo,"(' spin2:   ',3f15.6/' total:   ',3f15.6)")rmu1(2),rmu2(2),rhovxc(2),sum(rmu1),sum(rmu2),sum(rhovxc)
       write(stdo,"(' val*vef  ',3f15.6/' val chg: ',3f15.6)") vefv1,vefv2,valvef,qv1,qv2,qloc
       if(nsp == 2) write(stdo,"(' val mom: ',3f15.6,'    core:',f11.6)") a1,a2,aloc,alocc
       write(stdo,"(' core chg:',3f15.6)") qcor1,qcor2,qlocc
    endif
    call tcx('locpt2')
  end subroutine locpt2

  subroutine elfigr(nc,stdo,z,efg1)
    !use m_mathlib,only: htridi,imtql2,htribk
    !- Computation of electric field gradient
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nc    :number of classes or sites
    !i   stdo  :standard output
    !i   z     :nuclear charge
    ! o Inputs/Outputs
    ! o  efg1  :input:  l=2 part of electric field at nucleus
    ! o        :output: Electric field gradient
    !l Local variables
    !l         :
    !r Remarks
    !r
    !u Updates
    !u   30 Sep 08 Adapted from old FP (A Svane)
    ! ----------------------------------------------------------------------
    implicit none
    integer :: nc,stdo
    double precision :: z(nc),efg1(5,1)
    integer :: i,ifi,ic,ifesn,j,n,ier,imax,ixx,iyy
    double precision :: v(3,3),vi(3,3),d(3),e(3),e2(3),tau(2,3)
    double precision :: tr(3,3),ti(3,3),da(3)
    double precision :: conv1,emmsfe,emmssn,Qfe,Qsn,pi,f0,s3,conv2,conv3, &
         ax,bx,cx,dx,ex,dmax,eta,split
    data conv1 /162.1d0/
    data emmsfe,emmssn /4.8d-8,7.97d-8/
    data Qfe,Qsn /0.21d-28,-0.08d-28/
    pi = 4d0*datan(1d0)
    f0 = -dsqrt(15d0/16d0/pi)
    s3 = dsqrt(3d0)
    conv2 = conv1*3d0
    do  i = 1, 1
       ifi = stdo
       if(i == 2) ifi=71
       write(ifi,'(/," Electric Field Gradient:")')
       write(ifi,'(" Site ",5x,"Tensor axes",5x,"esu/cm^2",4x,"V/m^2 ",4x,"eta",4x,"line splitting")')
       write(ifi,'("      ",5x,"           ",5x," x10^13 ",4x,"x10^19",4x,"   ",4x,"    (mm/s)")')
    enddo
    do  ic = 1, nc
       if (z(ic) > 0.01d0) then
          ifesn = 1
          if (dabs(z(ic)-26d0) < 0.01d0) then
             conv3 = conv2*1.0d19*0.5d0*Qfe/emmsfe
          else if (dabs(z(ic)-50d0) < 0.01d0) then
             conv3 = conv2*1.0d19*0.5d0*Qsn/emmssn
          else
             ifesn = 0
             conv3 = 0d0
          endif
          ax = efg1(1,ic)
          bx = efg1(2,ic)
          cx = efg1(3,ic)
          dx = efg1(4,ic)
          ex = efg1(5,ic)
          v(1,1) =  2d0*ex-2d0/s3*cx
          v(2,2) = -2d0*ex-2d0/s3*cx
          v(3,3) = 4d0/s3*cx
          v(1,2) = 2d0*ax
          v(1,3) = 2d0*dx
          v(2,3) = 2d0*bx
          v(2,1) = v(1,2)
          v(3,1) = v(1,3)
          v(3,2) = v(2,3)
          do  i = 1, 3
             do  j = 1, 3
                v(i,j) = f0*v(i,j)
                vi(i,j) = 0d0
                tr(i,j) = 0d0
                ti(i,j) = 0d0
             enddo
             tr(i,i) = 1d0
          enddo
          n = 3
          diag: block
            complex(8):: tt(3,3),oo(3,3),vv(3,3),img=(0d0,1d0)
            integer:: nmx=3,nev
            oo(:,1) = [complex(8):: 1d0,0d0,0d0]
            oo(:,2) = [complex(8):: 0d0,1d0,0d0]
            oo(:,3) = [complex(8):: 0d0,0d0,1d0]
            vv=v+img*vi
            call diagcv(oo,vv,tt,n,d,nmx,1d99,nev)
            tr=dreal(tt)
          endblock diag
          !     write(stdo,99) ic
          ! 99  format(/' EFG, class = ',I3)
          do  i = 1, 3
             da(i) = dabs(d(i))
          enddo
          dmax = 0d0
          imax = 0
          do  i = 1, 3
             if (da(i) > dmax) then
                dmax = da(i)
                imax = i
             endif
          enddo
          !  d(imax) is field gradient (Vzz)
          ixx = mod(imax,3)+1
          iyy = mod(imax+1,3)+1
          eta = 0d0
          if(dabs(d(imax)) > 1.d-2) eta = dabs((d(ixx)-d(iyy))/d(imax))
          !     do  i = 1, 3
          !        write(stdo,98) conv1*d(i),conv2*d(i),(tr(j,i),j=1,3)
          !     enddo
          !     if (ifesn.eq.0) write(stdo,97) eta
          split = conv3*da(imax)*dsqrt(1d0+eta**2/3d0)
          !     if (ifesn.eq.1) write(stdo,96) eta,split
          ! 98  format(3X,F12.4,' 10**13 esu/cm**2',3X,F12.4,' 10**19 V/m*2',
          !    +         /,12X,3F12.6)
          ! 97  format(/' eta = ',F12.4/)
          ! 96  format(/' eta = ',F12.4,' line splitting = ',F12.4,' mm/s'/)
          do  i = 1, 3
             if(i == 1) then
                write(ifi,'(i4,3x,3f6.3,2x,f8.2,2x,f8.2,5x,f6.3,5x,f8.5)') &
                     ic,(tr(j,i),j=1,3),conv1*d(i),conv2*d(i),eta,split
             else
                write(ifi,'(4x,3x,3f6.3,2x,f8.2,2x,f8.2,5x)') &
                     (tr(j,i),j=1,3),conv1*d(i),conv2*d(i)
             endif
          enddo
          write(ifi,'(1x)')
       endif
    enddo
  end subroutine elfigr
  !$$$!!---------------------
  !$$$      subroutine adddipole(v,rofi,nr,nlml,nsp, idipole, basr)
  !$$$!! v=v+x (or y or z). to calculate dipole matrix of Hamiltonian
  !$$$      integer:: nr,nlml,nsp,idipole,isp,nrmx
  !$$$      real(8):: v(nr,nlml,nsp),rofi(nr),basr
  !$$$      real(8),parameter:: pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi),srfpi3 = dsqrt(4d0*pi/3d0)
  !$$$      do isp=1,nsp
  !$$$        v(1:nr,1,isp) = v(1:nr,1,isp) + basr * srfpi
  !$$$        if(idipole==1) v(1:nr,4,isp) = v(1:nr,4,isp)+rofi(1:nr)* srfpi3 !x
  !$$$        if(idipole==2) v(1:nr,2,isp) = v(1:nr,2,isp)+rofi(1:nr)* srfpi3 !y
  !$$$        if(idipole==3) v(1:nr,3,isp) = v(1:nr,3,isp)+rofi(1:nr)* srfpi3 !z
  !$$$      enddo
  !$$$      end
end module m_locpot

