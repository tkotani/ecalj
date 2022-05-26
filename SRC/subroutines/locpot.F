      module m_locpot
      use m_lgunit,only:stdo
      public locpot
      private
      contains
      subroutine locpot (       !lcplxp , 
     .      sv_p_orhoat , qmom , vval , gpot0 , sv_p_osig , sv_p_otau  !slat , 
     .    , sv_p_oppi, ohsozz,ohsopm, ppnl , hab , vab , sab , vvesat , cpnvsa , rhoexc 
     .    , rhoex , rhoec , rhovxc , rvepsv , rvexv , rvecv , rvvxcv , 
     .    rveps , rvvxc , valvef , xcore , focexc , focex , focec , focvxc 
     .    , sqloc , sqlocc , saloc , qval , qsc , job , rhobg , 
     .     nlibu , lmaxu , vorb , lldau,novxc)!,idipole )
      use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,nkaph,lxcf,lhh,nkapii,nkaphh,alat
      use m_lmfinit,only: n0,nppn,nab,nrmx,nkap0,nlmx,nbas,nsp,lso,ssite=>v_ssite, sspec=>v_sspec,mxcst4
      use m_MPItk,only: master_mpi
      use m_struc_def  
      use m_uspecb,only:uspecb
      use m_ftox
c      use m_iors,only: v1pot
C- Make the potential at the atomic sites and augmentation matrices.
C ----------------------------------------------------------------------
Ci Inputs
Ci   orhoat:vector of offsets containing site density
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   vval  :electrostatic potential at MT boundary; needed
Ci         :to computed matrix elements of local orbitals.
Ci   gpot0 :integrals of local gaussians * phi0~
Ci         :phi0~ is the estatic potential of the interstitial
Ci   job   :1s  digit
Ci         : 1 make core and augmentation matrices
Ci         :10s digit
Ci         : 0 update potential used to define basis,
Ci         :   provided 1s digit job is also set
Ci         : 1 do not update potential
Ci         :   NB: Caller can also suppress update for specific species
Ci         :   through 4's bit of species->mxcst.
Ci         :100s digit
Ci         : 1 exclude exchange-correlation potential in making
Ci         :   matrix elements of the potential.
Ci         :   Also implies floating potential.
Ci         :1000s digit
Ci         :1 Make rveps and rvvxc
Ci         :10000s digit
Ci         :1 write sphere density for site ib to file rhoMT.ib
Ci         :100000s digit
Ci         :1 write sphere density for site ib to file rhoMT.ib
Ci  rhobg  :compensating background density
Ci  lso    :if nonzero, calculate LzSz matrix elements
Ci  lcplxp=1 only now :0 if ppi is real; 1 if ppi is complex
Ci  nlibu  : max number of lda+u blocks
Ci  lmaxu  : max l for U
Ci  vorb   : orbital dependent potential matrices
Ci  lldau  :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat beginning at dmats(*,lldau(ib))
Co Outputs
Co   ... the following are summed over all spheres
Co   vvesat:integral of valence density times electrostatic potential
Co   cpnvsa:integral of core+nucleus times electrostatic potential
Co          = int [ (rhoc-z) ves1~ - (rhocsm-rhonsm) ves2~ ]
Co          where ves1~ is the estat potential of the true density
Co          where ves2~ is the estat potential of the smooth density
Co   rhoexc:integral of density times xc energy density
Co   rhoex :integral of density times exch. energy density
Co   rhoec :integral of density times corr. energy density
Co   rhovxc:integral of density times xc potential
Co   rvepsv:integral of valence density times exc(valence density)
Co   rvexv :integral of valence density times ex(valence density)
Co   rvecv :integral of valence density times ec(valence density)
Co   rvvxcv:integral of valence density times vxc(valence density)
Co   rveps :int rhov*exc(rhotot) (only made if 1000's digit job set)
Co   rvvxc :int rhov*vxc(rhotot) (only made if 1000's digit job set)
Co   valvef:integral (rho1*(v1-2Z/r) - rho2*vsm) ??
Co   xcore :integral rhoc*(v1-2Z/r)
Co   focexc:integral of smoothed core and xc energy
Co   focex :integral of smoothed core and exchange energy
Co   focec :integral of smoothed core and correlation energy
Co   focvxc:integral of smoothed core and xc potential
Co   sqloc :total valence charge rho1-rho2 in sphere
Co   saloc :total valence magnetic moment in sphere
Co         :+core moment for sites with core hole
Co   sqlocc:total core charge qcor1-qcor2 in sphere
Co   qval  :nominal total valence charge qv-z
Co   qsc   :semicore charge from local orbitals
Co   osig  :augmentation overlap integrals; see augmat
Co   otau  :augmentation kinetic energy integrals; see augmat
Co   oppi  :augmentation kinetic + potential integrals; see augmat
Co   ppnl  :NMTO-like potential parameters
Co   hab   :integrals of the ham. with true w.f.  See Remarks in augmat
Co   vab   :integrals of the pot. with true w.f.  See Remarks in augmat
Co   sab   :integrals of    unity with true w.f.  See Remarks in augmat
Cl Local variables
Cl   lfltwf:T  update potential used to define basis
Cl   iblu  :index to current LDA+U block
Ci   idu   :idu(l+1)=1 => this l has a nonlocal U matrix
Cr Remarks
Cu Updates
Cu   10 Sep 08 Added electric field gradient (AxSv)
Cu   02 Jan 06 adds core magnetic moment to saloc
Cu   09 Nov 05 Convert dmat to complex form
Cu   06 Jul 05 Parameters for extended local orbitals are now input
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   12 Jun 05 Optionally write sphere densities to file
Cu   27 Apr 05 Added LDA+U  (Lambrecht)
Cu    1 Sep 04 Adapted mkpot to handle complex ppi; fold so into ppi
Cu   15 Jul 04 First implementation of extended local orbitals
Cu   14 Mar 04 Makes rhov*exc and rhov*vxc
Cu   14 Mar 04 Makes rhov*exc and rhov*vxc
Cu   19 Sep 03 (ATP) Enabled partial core occupation
Cu   02 Oct 02 (WRL) Added background potential
Cu    9 May 02 Added species-specific freezing of potential
Cu    8 May 02 Added rhoex and rhoec (T. Miyake)
Cu    7 May 02 Added rvexv and rvecv (T. Miyake)
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   17 Sep 01 Returns qsc.  Altered argument list.
Cu   28 Aug 01 Extended to local orbitals.  Altered argument list.
Cu   15 Aug 01 Generates rvepsv and rvvxcv.  New argument list
Cu   18 Apr 01 Added ability to exclude exchange-correlation potential
Cu   20 Feb 01 Added ppnl to potential parameters generated
Cu   13 Jun 00 spin polarized
Cu    1 May 00 Adapted from nfp locpot.f
Cu   17 Jun 98 Adapted from DLN to run parallel on SGI
C ----------------------------------------------------------------------
      implicit none
      integer::  job , i_copy_size,ibx,ir,isp,l,lm !lcplxp,
      type(s_rv1) :: sv_p_orhoat(3,nbas)
      type(s_cv1) :: sv_p_oppi(3,nbas)
      type(s_sblock):: ohsozz(3,nbas),ohsopm(3,nbas)
      type(s_rv1) :: sv_p_otau(3,nbas)
      type(s_rv1) :: sv_p_osig(3,nbas)
      integer nlibu,lmaxu,lldau(nbas),iblu
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
c      parameter (n0=10,nppn=12,nab=9)
      real(8):: qmom(1) , vval(1)
c      type(s_site)::ssite(nbas)
c      type(s_spec)::sspec(*)
      double precision cpnvsa,rhoexc(2),rhoex(2),rhoec(2),rhovxc(2),
     .focexc(2),focex(2),focec(2),focvxc(2),qval,sqloc,sqlocc,saloc,
     .valvef,vvesat,rvepsv,rvexv,rvecv,rvvxcv,xcore,rveps,rvvxc,
     .hab(nab,n0,nsp,nbas),vab(nab,n0,nsp,nbas),sab(nab,n0,nsp,nbas),
     .gpot0(1),ppnl(nppn,n0,nsp,nbas),rhobg
C ... Local parameters
      character spid*8
c      integer nlmx
c      parameter (nlmx=64)
      integer lh(nkap0),nkapi,nkape,k
      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      double precision ehl(n0),rsml(n0)
      double precision pi,srfpi,y0,rofi(nrmx),rwgt(nrmx),pnu(n0,2),
     .pnz(n0,2),gpotb(81),z,a,rmt,rsma,qc,ceh,rfoc,
     .qcorg,qcorh,qsc,cofg,cofh,qsca,rg,qv,valvs,cpnvs,
     .rhexc(2),rhex(2),rhec(2),rhvxc(2),qloc,qlocc,valvt,xcor,
     .fcexc(2),fcex(2),fcec(2),fcvxc(2),aloc,alocc,rvepvl,rvexl,
     .rvecl,rvvxvl,rveptl,rvvxtl
C ... for sm. Hankel tails
      integer iltab(nbas)
      double precision rs3,vmtz
C ... for LDA+U
      integer idu(4)
C ... for core hole
      character chole*8
      integer kcor,lcor
      double precision qcor(2),qc0,qsc0
      real(8),allocatable:: wk(:),rhol1(:),rhol2(:),
     .v1(:),v2(:),v1es(:),v2es(:),efg(:,:),zz(:)
      logical lfltwf
      integer:: i , nglob , ipr , iprint , j1 , ib , is , lmxl 
     .     , lmxa , nr , lmxb , kmax , lfoc , nrml , nlml,ifivesint,ifi
      logical,save:: secondcall=.false.
      logical :: phispinsym,cmdopt0,readov0,v0fix,v0write,novxc
      integer::lsox
      character*20::strib
      character strn*120
      real(8):: ov0mean,pmean
c      integer:: idipole
      call tcn('locpot')
      ipr = iprint()
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0 = 1d0/srfpi
      k = nrmx*nlmx*nsp
      allocate(wk(k),rhol1(k),rhol2(k),v1(k),v2(k),v1es(k),v2es(k))
      allocate(efg(5,nbas),zz(nbas))
      if (ipr .ge. 30) write(stdo,"(/' locpot:')")
      call iinit(iltab,nbas)
      vvesat = 0d0
      cpnvsa = 0d0
      do  i = 1, 2
        rhoexc(i) = 0d0
        rhoex(i)  = 0d0
        rhoec(i)  = 0d0
        rhovxc(i) = 0d0
        focexc(i) = 0d0
        focex(i)  = 0d0
        focec(i)  = 0d0
        focvxc(i) = 0d0
      enddo
      rvepsv  = 0d0
      rvexv   = 0d0
      rvecv   = 0d0
      rvvxcv  = 0d0
      rveps   = 0d0
      rvvxc   = 0d0
      valvef  = 0d0
      xcore   = 0d0
      qval    = 0d0
      qsc     = 0d0
      sqloc   = 0d0
      sqlocc  = 0d0
      saloc   = 0d0
      j1 = 1
      iblu = 0
      ifivesint=-999
      if(master_mpi) open(newunit=ifivesint,file='vesintloc',form='formatted',status='unknown')
      do  ib = 1, nbas
        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
        z=sspec(is)%z
        qc=sspec(is)%qc
        rg=sspec(is)%rg
        spid=sspec(is)%name
        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        rsma=sspec(is)%rsma
        lmxa=sspec(is)%lmxa
        lmxl=sspec(is)%lmxl
        lmxb=sspec(is)%lmxb
        zz(ib)=z
        if (lmxa .eq. -1) cycle !goto 10 !floating orbital
        i_copy_size=size(sspec(is)%idu)
        call icopy(i_copy_size,sspec(is)%idu,1,idu,1)
        kmax=sspec(is)%kmxt
c        i=sspec(is)%mxcst
C       Float wave functions if:
C         100's digit job > 0 set   OR:
C         4's bit mxcst=0   AND   10's digit job=0  AND  1's digit job=1
c        lfltwf = mod(job,1000) .ge. 100 .or.
        lfltwf = !mod(job,1000) .ge. 100 .or. !aug2013 comment out for lmfgw with FRZWF=T
     .  (.not.mxcst4(is)) .and. mod(job/10,10).eq.0 .and. mod(job,10).eq.1
c     .  (.not.mxcst4(is)) mod(i/4,2).eq.0 .and. mod(job/10,10).eq.0 .and. mod(job,10).eq.1
        call corprm(sspec,is,qcorg,qcorh,qsca,cofg,cofh,ceh,lfoc,rfoc,z)
        chole=sspec(is)%coreh
        call gtpcor(sspec,is,kcor,lcor,qcor)
        call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc0,qv,qsc0)
        if (qsc0 .ne. qsca .or. qc .ne. qc0-qsc0) then
          if(iprint()>0)
     &    write(stdo,ftox)' is=',is,'qsc0=',ftof(qsc0),'qsca',ftof(qsca),'qc',ftof(qc),'qc0',ftof(qc0)
!          =%d',' ',120,stdo,is,qsc0,qsca,qc,qc0)
          call rxs('problem in locpot -- possibly low LMXA or orbital mismatch, species ',spid)
        endif
        qval = qval+qv+qsca
        qsc  = qsc+qsca
        nlml = (lmxl+1)**2
        nrml = nr*nlml
        if (ipr >= 20) then
          write(stdo,200) ib,z,rmt,nr,a,nlml,rg,lfltwf
  200     format(/' site',i3,'  z=',f5.1,'  rmt=',f8.5,'  nr=',i3,
     .    '   a=',f5.3,'  nlml=',i2,'  rg=',f5.3,'  Vfloat=',l1)
          if (kcor .ne. 0) then
            if (qcor(1) .ne. 0 .or. qcor(2) .ne. 0) then
              call info5(30,0,0,' core hole:  kcor=%i  lcor=%i  qcor=%d  amom=%d',kcor,lcor,qcor,qcor(2),0)
            endif
          endif
        endif
        call rxx(nr .gt. nrmx,  'locpot: increase nrmx')
        call rxx(nlml .gt. nlmx,'locpot: increase nlmx')
        call radmsh(rmt,a,nr,rofi)
        call radwgt(rmt,a,nr,rwgt)
C     ... Write true density to file rhoMT.ib
        if (mod(job/10000,10) .eq. 1 .or. mod(job/10000,10) .eq. 3) then
          call wrhomt('rhoMT.','density',ib,sv_p_orhoat(1,ib)%v,rofi,nr,nlml,nsp)
        endif
C   --- Make potential and energy terms at this site ---
        call locpt2 (job , z , rmt , rg , a , nr , nsp , cofg , cofh 
     .   , ceh , rfoc , lfoc , nlml , qmom ( j1 ) , vval ( j1 ) , rofi 
     .   , rwgt , sv_p_orhoat( 1 , ib )%v , sv_p_orhoat( 2 , ib )%v , 
     .   sv_p_orhoat( 3 , ib )%v , rhol1 , rhol2 , v1 , v2 , v1es , v2es 
     .   , wk , valvs , cpnvs , rhexc , rhex , rhec , rhvxc , rvepvl , 
     .   rvexl , rvecl , rvvxvl , rveptl , rvvxtl , valvt , xcor , qloc 
     .   , qlocc , aloc , alocc , gpotb , fcexc , fcex , fcec , fcvxc 
     .       , rhobg , efg ( 1 , ib ),ifivesint) !,novxc )
!! write density 1st(true) component and counter components.
        if(cmdopt0('--density').and.master_mpi.and.secondcall) then
           write(stdo,"(' TotalValenceChange diff in MT;  ib,\int(rho2-rho1)=',i5,f13.5)") ib,qloc
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
C   ... Write true potential to file vtrue.ib
        if (mod(job/10000,10) .eq. 2 .or. mod(job/10000,10) .eq. 3) then
          call wrhomt('vtrue.','potential',ib,v1,rofi,nr,nlml,nsp)
        endif
C   ... Update the potential used to define basis set
        if (lfltwf) then
          do i = 0, nsp-1
            call dpscop(v1, ssite(ib)%rv_a_ov0, nr, 1+nr*nlml*i, 1+nr*i, y0)
          enddo
        endif
!! spin averaged oV0 to generate phi and phidot. takaoAug2019
        phispinsym= cmdopt0('--phispinsym')
        if(phispinsym) then
           if(master_mpi.and.nsp==2) then
             write(6,*) 'locpot: --phispinsym mode: use spin-averaged potential for phi and phidot'
           endif   
           do ir=1,nr
              ov0mean = 0d0
              do isp=1,nsp
                 ov0mean = ov0mean + ssite(ib)%rv_a_ov0( ir + nr*(isp-1) )
              enddo
              ov0mean = ov0mean/nsp
              do isp=1,nsp
                 ssite(ib)%rv_a_ov0(ir + nr*(isp-1))= ov0mean
              enddo
           enddo
        endif
           
ccccccccccccccccccc experimental case --v0fix
        v0fix= cmdopt0('--v0fix')
        if(v0fix) then
           inquire(file='v0pot.'//char(48+ib),exist=readov0)
           if(readov0) then
              block 
              real(8):: ov0(nr)
              open(newunit=ifi,file='v0pot.'//char(48+ib),form='unformatted')
              read(ifi) ov0(1:nr)
              close(ifi)
              do ir=1,nr
              do isp=1,nsp
                 ssite(ib)%rv_a_ov0(ir + nr*(isp-1))= ov0(ir)
              enddo
              enddo
              end block
           else
              call rx('no v0pot files')
           endif
        endif   
cccccccc
        v0write= cmdopt0('--v0write')
        if(v0write) then
        do ir=1,nr
           ov0mean = 0d0
           do isp=1,nsp
              ov0mean = ov0mean + ssite(ib)%rv_a_ov0( ir + nr*(isp-1) )
           enddo
           ov0mean = ov0mean/nsp !spin averaged
           do isp=1,nsp
              ssite(ib)%rv_a_ov0(ir + nr*(isp-1))= ov0mean
           enddo
        enddo
        open(newunit=ifi,file='v0pot.'//char(48+ib),form='unformatted')
        write(ifi) ssite(ib)%rv_a_ov0(1:nr)
        close(ifi)
        if(ib==nbas) readov0=.true.
        endif      
ccccccccccccccccccccccccccccccccccccccc
        
        if(master_mpi.and.nsp==2)then
           do l=0,lmxa
              write(6,"('  ibas l=',2i3,' pnu(1:nsp) pnz(1:nsp)=',4f10.5)") ib,l,pnu(l+1,1:nsp),pnz(l+1,1:nsp)
           enddo
        endif
C   ... Store the potential used in mkrout to calculate the core
        do  i = 0, nsp-1
           call dpscop ( v1 , ssite(ib)%rv_a_ov1, nr , 1 + nr * nlml * i , 1 + nr* i , y0 )
c           v1pot(ib)%v=y0*v1
        enddo
C   ... Accumulate terms for LDA total energy
        vvesat = vvesat + valvs
        cpnvsa = cpnvsa + cpnvs
        valvef = valvef + valvt
        rvepsv = rvepsv + rvepvl
        rvexv  = rvexv  + rvexl
        rvecv  = rvecv  + rvecl
        rvvxcv = rvvxcv + rvvxvl
        rveps  = rveps  + rveptl
        rvvxc  = rvvxc  + rvvxtl
        sqloc  = sqloc  + qloc
        saloc  = saloc  + aloc
        sqlocc = sqlocc + qlocc
        rhoexc(1) = rhoexc(1) + rhexc(1)
        rhoex(1)  = rhoex(1)  + rhex(1)
        rhoec(1)  = rhoec(1)  + rhec(1)
        rhovxc(1) = rhovxc(1) + rhvxc(1)
        focexc(1) = focexc(1) + fcexc(1)
        focex(1)  = focex(1)  + fcex(1)
        focec(1)  = focec(1)  + fcec(1)
        focvxc(1) = focvxc(1) + fcvxc(1)
        if (nsp .eq. 2) then
          rhoexc(2) = rhoexc(2) + rhexc(2)
          rhoex(2)  = rhoex(2)  + rhex(2)
          rhoec(2)  = rhoec(2)  + rhec(2)
          rhovxc(2) = rhovxc(2) + rhvxc(2)
          focexc(2) = focexc(2) + fcexc(2)
          focex(2)  = focex(2)  + fcex(2)
          focec(2)  = focec(2)  + fcec(2)
          focvxc(2) = focvxc(2) + fcvxc(2)
        endif
        if (lfoc==0) xcore = xcore + xcor
C       Check for core moment mismatch ; add to total moment
        if (kcor/=0) then
          if (dabs(qcor(2)-alocc) .gt. 0.01d0) then
            call info5(10,0,0,' (warning) core moment mismatch spec %i:'
     .      //'  input file=%;6d  core density=%;6;4d',
     .      is,qcor(2),alocc,0,0)
          endif
          saloc = saloc + qcor(2)
        endif
C   --- Make augmentation matrices sig, tau, ppi ---
        if (mod(job,10)==1) then
C     ... Smooth Hankel tails for local orbitals
          rsmh= 0d0 
          eh  = 0d0
          call uspecb(is,rsmh,eh)
          nkapi=nkapii(is)
          nkape=nkaphh(is)
          if (nkape > nkapi) then
            rsml=rsmh(:,nkaph) !call dcopy(n0,rsmh(1,nkaph),1,rsml,1)
            ehl =  eh(:,nkaph) !call dcopy(n0,eh(1,nkaph),  1,ehl, 1)
          endif
C     ... Use effective potentials with modified xc
          if (novxc) then 
            v1=v1es !call dcopy(nrml*nsp,v1es,1,v1,1)
            v2=v2es !call dcopy(nrml*nsp,v2es,1,v2,1)
          endif
!!
          lsox=lso
c          if(idipole/=0) then
c             call adddipole(v1,rofi,nr,nlml,nsp,idipole,ssite(ib)%pos(idipole)*alat)
c             call adddipole(v2,rofi,nr,nlml,nsp,idipole,ssite(ib)%pos(idipole)*alat)
c          endif
          if(.not.novxc .and. cmdopt0('--socmatrix') ) then
             lsox=1
          endif   
          
          if (ipr >= 20) write(stdo,467) y0*(gpot0(j1)-gpotb(1))
  467     format(' potential shift to crystal energy zero:',f12.6)
          call augmat ( z , rmt , rsma , lmxa , pnu , pnz , kmax , nlml
     .           , a , nr , nsp , lsox , rofi , rwgt , rv_a_ocg, iv_a_ojcg, iv_a_oidxcg
     .           , ssite(ib)%rv_a_ov0 , v1 , v2 , gpotb , gpot0 ( j1 ) , nkaph , nkapi , 
     .           lmxb , lhh(:,is) , eh , rsmh, ehl , rsml , rs3 , vmtz ,  lmaxu !lcplxp ,
     .         , vorb, lldau(ib), iblu, idu, sv_p_osig(1,ib) , sv_p_otau(1,ib),
     &           sv_p_oppi(1,ib),ohsozz(1,ib),ohsopm(1,ib), ppnl(1,1,1,ib) , 
     .           hab(1,1,1,ib),vab (1,1,1,ib), sab(1,1,1,ib) )
        endif
        j1 = j1+nlml
c   10   continue
      enddo
      if(cmdopt0('--density').and.master_mpi) secondcall=.true.
      if(master_mpi) close(ifivesint)
C     Electric field gradient
      if (ipr .gt. 40) call elfigr(nbas,stdo,zz,efg)
      deallocate(efg,zz)
      deallocate(wk,rhol1,rhol2,v1,v2,v1es,v2es)
      call tcx('locpot')
      end subroutine locpot

c$$$!!---------------------
c$$$      subroutine adddipole(v,rofi,nr,nlml,nsp, idipole, basr)
c$$$!! v=v+x (or y or z). to calculate dipole matrix of Hamiltonian      
c$$$      integer:: nr,nlml,nsp,idipole,isp,nrmx
c$$$      real(8):: v(nr,nlml,nsp),rofi(nr),basr
c$$$      real(8),parameter:: pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi),srfpi3 = dsqrt(4d0*pi/3d0)
c$$$      do isp=1,nsp
c$$$        v(1:nr,1,isp) = v(1:nr,1,isp) + basr * srfpi
c$$$        if(idipole==1) v(1:nr,4,isp) = v(1:nr,4,isp)+rofi(1:nr)* srfpi3 !x
c$$$        if(idipole==2) v(1:nr,2,isp) = v(1:nr,2,isp)+rofi(1:nr)* srfpi3 !y
c$$$        if(idipole==3) v(1:nr,3,isp) = v(1:nr,3,isp)+rofi(1:nr)* srfpi3 !z
c$$$      enddo
c$$$      end
      
      subroutine locpt2(job,z,rmt,rg,a,nr,nsp,cofg,cofh,ceh,rfoc,lfoc,
     .nlml,qmom,vval,rofi,rwgt,rho1,rho2,rhoc,rhol1,rhol2,v1,v2,v1es,
     .v2es,wk,vvesat,cpnves,rhoexc,rhoex,rhoec,rhovxc,rvepsv,
     .rvexv,rvecv,rvvxcv,rveps,rvvxc,valvef,xcore,qloc,
     .qlocc,aloc,alocc,gpotb,focexc,focex,focec,focvxc,rhobg,efg,ifivesint)!,novxc)
      use m_lmfinit,only: lxcf,nlmx,nrmx
      use m_ftox
C- Makes the potential at one site, and associated energy terms.
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :1s digit - 100s digit not used here.
Ci         :1000s digit
Ci         :1 Make rveps and rvvxc
Ci   z     :nuclear charge
Ci   rmt   :augmentation radius
Ci   rg    :smoothing radius for compensating gaussians used to
Ci         :correct the multipole moments of local smooth density
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
Ci         :cofg = Y0 * qcorg
Ci   cofh  :coefficient to Hankel part of pseudocore density (corprm)
Ci   ceh   :energy of hankel function to fit core tail
Ci   rfoc  :smoothing radius for hankel head fitted to core tail
Ci   lfoc  :switch specifying treatment of core density.
Ci         :0 => core w.f. have val,slo = 0 at rmt
Ci         :1 => core included explicitly with valence
Ci         :2 => core included perturbatively
Ci   nlml  :L-cutoff for rho1,rho2
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   vval  :boundary condition for estat potential at MT boundary.
Ci         :See remarks.
Ci   rofi  :radial mesh points for density and potential
Ci   rwgt  :radial mesh weights to integrate function on radial mesh
Ci         :Integral f(r) dr = sum_i f_i * wt_i
Ci   rho1  :local valence density = sum_ilm rho1(ilm) * r**2  Y_L(ilm),
Ci         :on radial mesh.
Ci   rho2  :local smoothed valence density, defined as rho1
Ci         :Local atomic valence density is rho1-rho2
Ci   rhoc  :core density times 4*pi*r*r
Ci   wk    :work array of dimension nr*nlml*nsp
Ci   rhobg: compensating background density
Co Outputs
Co   cpnves:integral of core+nucleus times electrostatic potential
Co          = int [ (rhoc-z) ves1~ - (rhocsm-rhonsm) ves2~ ]
Co          where ves1~ is the estat potential of the true density
Co          where ves2~ is the estat potential of the smooth density
Co   vvesat:integral of valence density times electrostatic potential
Co   rhoexc:integral of valence density times xc energy density
Co   rhoex :integral of valence density times exch. energy density
Co   rhoec :integral of valence density times corr. energy density
Co   rhovxc:integral of valence density times xc potential
Co   rvepsv:integral of valence density exc(valence density)
Co   rvexv :integral of valence density ex(valence density)
Co   rvecv :integral of valence density ec(valence density)
Co   rvvxcv:integral of valence density vxc(valence density)
Co   valvef:integral valence density times full potential:
Co         :used to compute double-counting for kinetic energy.
Co         : = vefv1 - vefv2; defined in Local variables below
Co   xcore :integral rhoc*(v1-2Z/r)
Co   qloc  :total valence charge in sphere rho1 - rho2
Co   aloc  :total valence magnetic moment in sphere
Co   qlocc :total core charge in sphere qcor1 - qcor2
Co   alocc :total core magnetic moment in sphere
Co   gpotb :integrals of gaussians(radius rg) times smooth ves.
Co         :Here smooth ves is the electrostatic potential of the
Co         :compensated sm density rhol2, w/ b.c. ves(rmax)=0
Co         :This is a local analog of gpot0 generated for the smooth
Co         :density.  Here gpotb is generated in a Y_lm expansion
Co         :of the gaussian orbitals
Co   focexc:(lfoc=2): integral rhochs * vxc
Co         :otherwise, 0
Co   focex :integral of smoothed core and exchange energy
Co   focec :integral of smoothed core and correlation energy
Co   focvxc:If lfoc=2, integral rhochs * (dvxc/drho2 * rho2)
Co         :Otherwise, 0
Co   rhol1 :full true electron density, rho1 + rhoc
Co   rhol2 :full smooth density, i.e. uncompensated rho2 plus
Co         :compensating gaussians + pseudocore charge
Co   v1    :true potential at this site, excluding nuclear contribution
Co         :but including background contribution; see Remarks
Co         :See Remarks concerning boundary condition at rmt
Co   v2    :Total smooth potential including background potential
Co         : = ves[rhol2] + ...
Co         :   ... lfoc=0 : vxc(rho2)
Co         :       lfoc=1 : vxc(rho2+sm-rhoc)
Co         :       lfoc=2 : vxc(rho2) + dvxc/drho * sm-rhoc
Co         :Apart from differences in l-truncation in the XC potential
Co         :(i.e. true and local sm. densities are different),
Co         :v1(nr,1,1)-2*z/y0/rmt = v2(nr,1,1)
Co         :See Remarks concerning boundary condition at rmt
Co  efg    :l=2 potential at nucleus
Cl Local variables
Cl   rhol1 :full electron density, rho1 + rhoc
Cl   rhol2 :full smooth compensated density rho2 + gval + gcor + gnuc
Cl   rhochs:Part of smooth core density contains Hankel tails
Cl         :Used to make vxc[n2] and exc[n2]
Cl         :rhochs = srfpi*cofh*xi(0)*r*r
Cl   vefv1 :int (valence true density)  * (total true potential):
Cl         : = int n1 * v1~
Cl   vefv2 :int (valence smooth density)* (total smooth potential):
Cl         : = int (n2*v2~) +  sum_L qmom_L gpotb_L
Cl         :   [- perturbation (focvxc(1)+focvxc(2)) if lfoc=2]
Cl   qcor1 :true core charge inside rmt
Cl   qcor2 :smoothed core charge inside rmt
Cl   rhonsm:nuclear density smoothed into gaussian of width rg
Cl   rhocsm:core density smoothed into gaussian of width rg
Cr Remarks
Cr   On boundary conditions for the estat potential vval at the MT
Cr   boundary.  In the original formulation (see Springer book) the
Cr   the total hamiltonian and energy did not depend on the choice of
Cr   vval, as the contribution from true and smooth potentials exactly
Cr   cancel.  However, this is no longer true when local orbitals are
Cr   present that have no corresponding smooth part.
Cr
Cr   The potential due to uniform background density rhobg is added
Cr                 vbg= -(4pi/3)*(rmt^2-r^2)*rhobg to v1 and v2
Cr   in spherical components .
Cu Updates
Cu   30 Sep 08 Return l=2 potential at nucleus in efg (A Svane)
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   19 Sep 02 added uniform bkg potential interactions (WRL)
Cu    8 Feb 02 rhoex and rhoec (T. Miyake)
Cu   14 Jan 02 rvexv and rvecv (T. Miyake)
Cu   15 Aug 01 Generates rvepsv and rvvxcv
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,lfoc,nlml,job
      double precision z,rmt,rg,a,cofg,cofh,ceh,rfoc,xcore,qloc,qlocc,
     .aloc,alocc,rhoexc(2),rhoex(2),rhoec(2),rhovxc(2),focexc(2),
     .focex(2),focec(2),focvxc(2),valvef,vvesat,
     .cpnves,rvepsv,rvexv,rvecv,rvvxcv,rveps,rvvxc,rhobg
      double precision rofi(nr),rwgt(nr),
     .qmom(nlml),vval(nlml),gpotb(nlml),
     .rho1(nr,nlml,nsp),rhol1(nr,nlml,nsp),
     .rho2(nr,nlml,nsp),rhol2(nr,nlml,nsp),
     .v1(nr,nlml,nsp),v1es(nr,nlml,nsp),
     .v2(nr,nlml,nsp),v2es(nr,nlml,nsp),
     .wk(nr,nlml,nsp),rhoc(nr,nsp)
      double precision efg(5)
      integer ipr,iprint,ll,i,isp,ilm,l,lxcfun,nglob,nrml,
     .isw,isw2
      double precision rhochs(nrmx*2),rhonsm(nrmx),df(0:20),cof(nlmx),
     .rhocsm(nrmx),xi(0:40),tmp(2)
      double precision afoc,ag,b,cof0,fac,gnu,pi,qv1,qv2,qcor1,qcor2,
     .r,rep1(2),rep2(2),rep1x(2),rep2x(2),rep1c(2),rep2c(2),
     .rhves1,rhves2,rmu1(2),rmu2(2),rvs1,rvs2,
     .rvsm(2),rvtr(2),samh,sfac,sgpotb,srfpi,sum1,sum2,sumg,sumh,top,
     .ves1,vales1,vales2,vcpn1,vcpn2,vefc1,vefv1,vefv2,vesc1,vesc2,
     .vesn1,vesn2,vnucl,vsum,vtr,y0,ddot,a1,a2,smrhoc
      double precision qs(2)
      real(8), allocatable :: fl(:,:,:)
      logical:: debug=.false.
      real(8):: ves1int,ves2int
      integer:: ifivesint
c      logical:: novxc
      call tcn('locpt2')
c      stdo = lgunit(1)
      ipr = iprint()
      isw = mod(job/1000,10)
      isw2 = mod(job/100000,10)*10
      rveps = 0
      rvvxc = 0
      alocc = 0
      if (isw .ne. 0) then
        allocate (fl(nr,nlml,nsp+1))
      else
        allocate (fl(1,1,1))
      endif
      if (nr .gt. nrmx) call rxi('locpt2: need nrmx at least',nr)
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4*pi)
      y0 = 1d0/srfpi
      if (nlml .gt. nlmx) call rxi('locpt2: increase nlmx, need',nlml)
      call stdfac(20,df)
      b = rmt/(dexp(a*nr-a)-1)
      ag = 1d0/rg
      afoc = 1d0/rfoc
      nrml = nr*nlml
C --- Make core and nucleus pseudodensities ---
      fac = 4*pi*(ag*ag/pi)**1.5d0
C ... Renormalize gaussian
      sumg = 0d0
      do  i = 2, nr
        r = rofi(i)
        gnu = fac* r*r * dexp(-ag*ag*r*r)
        sumg = sumg + rwgt(i)*gnu
      enddo
      if (dabs(sumg-1d0) .gt. 1d-4) write(stdo,ftox)' locpot (warning): '
     .//'large gaussian, integral=',ftod(sumg)
      sfac = 1d0/sumg
      fac = fac*sfac
      sumh = 0d0
C     Smooth nuc. and core rho, sm Hankel portion, true & smooth core q
      do  i = 1, nr
        r = rofi(i)
        gnu = fac * r*r * dexp(-ag*ag*r*r)
        call hansmr(r,ceh,afoc,xi,1)
        rhonsm(i) = -z*gnu
        smrhoc = srfpi*cofh*xi(0)*r*r
        rhocsm(i) = srfpi*cofg*gnu + smrhoc
        rhochs(i) = smrhoc
        sumh  = sumh + rwgt(i)*rhochs(i)
C       qcor1 = qcor1 + rwgt(i)*rhoc(i)
C       qcor2 = qcor2 + rwgt(i)*rhocsm(i)
C       sum1 = sum1 + rwgt(i)*rhonsm(i)
      enddo
      samh = -y0*cofh*4d0*pi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
      if (ipr.ge.20 .and. dabs(samh) .gt. 1d-6)
     .write(stdo,ftox)' sm core charge=',ftof(sumh),'=(sphere)+(spillout)=',ftof(samh-sumh),ftof(samh)
c     .//' = %;6d',' ',80,stdo,sumh,samh-sumh,samh)
C --- rhol1 = full electron density = rho1 + rhoc ---
      if(debug) print *,'locpt2: 1111'
      do  isp = 1, nsp
        rhol1(:,:,isp)= rho1(:,:,isp)
        rhol1(:,1,isp)= rhol1(:,1,isp)+ y0*rhoc(:,isp)
      enddo
C --- rhol2 = full smooth compensated density rho2+gval+gcor+gnuc ---
C     gval : qmom * compensating gaussians
C     Distribute core+nuclear charge equally over spins
      if(debug) print *,'locpt2: 2222'
      do  isp = 1, nsp
        do  ilm = 1, nlml
          l = ll(ilm)
          cof(ilm) = qmom(ilm)*4*pi/df(2*l+1)
          fac = sfac*(ag*ag/pi)**1.5d0 * (2*ag*ag)**l
          do  i = 1, nr
            r = rofi(i)
            gnu = cof(ilm)*fac* r**(l+2) * dexp(-ag*ag*r*r)
            rhol2(i,ilm,isp) = rho2(i,ilm,isp) + gnu/nsp
          enddo
        enddo
        do  i = 1, nr
          rhol2(i,1,isp) = rhol2(i,1,isp) + y0/nsp*(rhocsm(i)+rhonsm(i))
        enddo
      enddo
C ... Combine separate spin densities for electrostatics
      call splrho(0,nsp,nr,nlml,rho1,rho2,rhoc)
      call splrho(20,nsp,nr,nlml,rhol1,rhol2,rhoc)
C ... Add background density to spherical rhol1 and rhol2
      if(debug) print *,'locpt2: 3333'
      do i = 1, nr
        rhol1(i,1,1) = rhol1(i,1,1) + srfpi*rhobg*rofi(i)**2
        rhol2(i,1,1) = rhol2(i,1,1) + srfpi*rhobg*rofi(i)**2
      enddo
C ... Sphere charges; also check sphere neutrality for safety
      qv1   = srfpi*ddot(nr,rwgt,1,rho1,1)
      qv2   = srfpi*ddot(nr,rwgt,1,rho2,1)
      a1    = srfpi*ddot(nr,rwgt,1,rho1(1,1,nsp),1)
      a2    = srfpi*ddot(nr,rwgt,1,rho2(1,1,nsp),1)
      qcor1 =       ddot(nr,rwgt,1,rhoc,1)
      qcor2 =       ddot(nr,rwgt,1,rhocsm,1)
      qlocc = qcor1-qcor2
      if (nsp .eq. 2) then
        alocc =        ddot(nr,rwgt,1,rhoc(1,2),1)
      endif
      qloc  = qv1-qv2
      sum1  = srfpi*ddot(nr,rwgt,1,rhol1,1) - z
      sum2  = srfpi*ddot(nr,rwgt,1,rhol2,1)
      qlocc = qcor1-qcor2
      qloc  = qv1-qv2
      aloc  = a1-a2
      if (nsp .eq. 1) aloc = 0
      if (dabs(sum1-sum2)>1d-6)call rx1('locpt2: sphere not neutral: charge = %d',sum1-sum2)
      
c      if(.not.novxc) then
C --- Solve Poisson equation for the true and smooth densities ---
C     call dpzero(vval, nlml)
C     v1 = Ves[rho1]: true ES pot without nuclear contribution
         if(debug) print *,'locpt2: 4444'
         call poinsp(z,vval,nlml,a,b,v1,rofi,rhol1,wk,nr,rvs1,rhves1,  vnucl,vsum)
         if (nlml.ge.9.and.z.gt.0.01) then
            do ilm=5,9
               efg(ilm-4)=v1(5,ilm,1)/rofi(5)**2
            enddo
         else
            do ilm=1,5
               efg(ilm)=0.d0
            enddo
         endif
C     v2 = Ves[rhol2 = rho2+gval+gcor+gnuc] ---
         call poinsp(0d0,vval,nlml,a,b,v2,rofi,rhol2,wk,nr,rvs2,rhves2, vnucl,vsum)
c      else
c         v1=0d0
c         v2=0d0
c         rvs1=0d0
c         rvs2=0d0
c         rhves1=0d0
c         rhves2=0d0
c         vnucl=0d0
c         vsum=0d0
c      endif
      
C --- gpotb = integrals of compensating gaussians times smooth ves ---
      if(debug) print *,'locpt2: 5555'
      sgpotb = 0d0
      do  ilm = 1, nlml
        l = ll(ilm)
        cof0 = 4*pi/df(2*l+1)
        fac = sfac*(ag*ag/pi)**1.5d0 * (2*ag*ag)**l
        sum1 = 0d0
        do  i = 1, nr
          r = rofi(i)
          gnu = cof0*fac* r**(l+2) * dexp(-ag*ag*r*r)
          sum1 = sum1 + rwgt(i)*v2(i,ilm,1)*gnu
        enddo
        gpotb(ilm) = sum1
        sgpotb = sgpotb + qmom(ilm)*gpotb(ilm)
      enddo
C --- Electrostatic integrals involving spherical terms only ---
      vesn2 = 0d0
      vesc2 = 0d0
      vesc1 = 0d0
      vnucl = 0d0
      ves1int=0d0
      ves2int=0d0
      do  i = 2, nr
        ves1 = y0*v1(i,1,1) - 2*z/rofi(i)
        vesc1 = vesc1 + rwgt(i)*rhoc(i,1)*ves1
        vesn2 = vesn2 + rwgt(i)*rhonsm(i)*y0*v2(i,1,1)
        vesc2 = vesc2 + rwgt(i)*rhocsm(i)*y0*v2(i,1,1)
        vnucl = vnucl + rwgt(i)*rhol1(i,1,1)*(1d0/rofi(i)-1d0/rmt)
        ves1int = ves1int + rwgt(i)*y0*v1(i,1,1)*rofi(i)**2
        ves2int = ves2int + rwgt(i)*y0*v2(i,1,1)*rofi(i)**2
      enddo
      ves1int = 4d0*pi*(ves1int - z*rofi(nr)**2)
      ves2int = 4d0*pi*ves2int 
      if(ifivesint>0) then
      write(ifivesint,"(3f23.15,a)")ves1int-ves2int,ves1int,ves2int,' ! vesint1-vesint2 ves1int ves2int'
      endif
      vnucl = 2*srfpi*vnucl + 2*z/rmt + y0*vval(1)
      vesn1 = -z*vnucl
C ... Valence density times electrostatic potential
      vales1 = rvs1-vesc1
      vales2 = rvs2-vesn2-vesc2
      vvesat = vales1-vales2
C ... Core plus nucleus times estatic potential
      vcpn1  = vesc1 + vesn1
      vcpn2  = vesn2 + vesc2
      cpnves = vcpn1 - vcpn2
C ... Subtract background before doing exchange correlation
      do i=1,nr
        rhol1(i,1,1)=rhol1(i,1,1)-srfpi*rhobg*rofi(i)**2
        rhol2(i,1,1)=rhol2(i,1,1)-srfpi*rhobg*rofi(i)**2
      enddo
C ... Restore separate spin densities; copy estat pot to spin2
      if (nsp .eq. 2) then
        call splrho(1,nsp,nr,nlml,rho1,rho2,rhoc)
        call splrho(21,nsp,nr,nlml,rhol1,rhol2,rhoc)
        call dcopy(nrml,v1,1,v1(1,1,2),1)
        call dcopy(nrml,v2,1,v2(1,1,2),1)
      endif
C ... Preserve potentials without exchange-correlation
      call dcopy(nrml*nsp,v1,1,v1es,1)
      call dcopy(nrml*nsp,v2,1,v2es,1)
C ... Generate valence-only rvepsv and rvvxcv (uses v1 as work array)
      call pshpr(max(ipr-31,min(ipr,10)))
      lxcfun = lxcf
      if(debug) print *,'locpt2: 55551111 isw2 lxcfun=',isw2,lxcfun
      if(debug) write(6,'(a)')' === rho1 valence true density ==='
      call vxcnsp(isw2,a,rofi,nr,rwgt,nlml,nsp,rho1,lxcfun,0,0,0,0,0,
     .rep1,rep1x,rep1c,rmu1,v1,fl,qs)
      if(debug) print *,'locpt2: 55551111aaaa'
      if(debug) write(6,'(a)')' === rho2 valence counter density ==='
      call vxcnsp(isw2,a,rofi,nr,rwgt,nlml,nsp,rho2,lxcfun,0,0,0,0,0,
     .rep2,rep2x,rep2c,rmu2,v1,fl,qs)
      rvvxcv = rmu1(1) - rmu2(1) + rmu1(2) - rmu2(2)
      rvepsv = rep1(1) - rep2(1) + rep1(2) - rep2(2)
      rvexv  = rep1x(1) - rep2x(1) + rep1x(2) - rep2x(2)
      rvecv  = rep1c(1) - rep2c(1) + rep1c(2) - rep2c(2)
      call dcopy(nrml*nsp,v1es,1,v1,1)
      call poppr
      if(debug) print *,'locpt2: 5555222'
C --- Add xc potentials to v1 and v2 ---
      call pshpr(max(ipr-11,min(ipr,10)))
      lxcfun = lxcf !globalvariables%lxcf
      call dpzero(focexc,2)
      call dpzero(focex,2)
      call dpzero(focec,2)
      call dpzero(focvxc,2)
      call info0(30,0,0,' Exchange for true density:')
      if(debug) write(6,'(a)')' === rhol1 valence+core density ==='
      call vxcnsp(isw+isw2,a,rofi,nr,rwgt,nlml,nsp,rhol1,lxcfun,
     .0,0,0,0,0,rep1,rep1x,rep1c,rmu1,v1,fl,qs)
      if (isw .eq. 1) then
        call dpzero(tmp,2)
        call vxcns5(1,31,'rhov*vxc',nlml,nsp,nr,rofi,rwgt,rho1,fl,xi,tmp)
        rvvxc = tmp(1) + tmp(2)
        do  isp = 1, nsp
          call vxcns5(1,31,'rhov*exc',nlml,1,nr,rofi,rwgt,rho1(1,1,isp),
     .    fl(1,1,3),xi,tmp(isp))
        enddo
        rveps = tmp(1) + tmp(2)
      endif
      if(debug) print *,'locpt2: 5555333'

C     call prrmsh('v1',rofi,v1,nr,nr,nlml)
C     If no core treatment v2 += vxc(rho2)
      if(debug)write(6,'(a)')' === rho2 ->valence+smooth core density ==='
      if (lfoc .eq. 0) then
        call vxcnsp(isw+isw2,a,rofi,nr,rwgt,nlml,nsp,rho2,lxcfun,
     .  0,0,0,0,0,rep2,rep2x,rep2c,rmu2,v2,fl,qs)
C     Otherwise v2 += vxc(rho2 + sm core), directly or perturbatively:
      else if (lfoc .eq. 1) then
        if (ipr .gt. 40) print *, 'exchange for smooth density, foca=1:'
        do  isp = 1, nsp
c          call dpadd(rho2(1,1,isp),rhochs,1,nr,y0/nsp)
          rho2(1:nr,1,isp)=rho2(1:nr,1,isp) + y0/nsp*rhochs(1:nr)
        enddo
        call vxcnsp(isw+isw2,a,rofi,nr,rwgt,nlml,nsp,rho2,lxcfun,
     .  0,0,0,0,0,rep2,rep2x,rep2c,rmu2,v2,fl,qs)
        do  isp = 1, nsp
c          call dpadd(rho2(1,1,isp),rhochs,1,nr,-y0/nsp)
          rho2(1:nr,1,isp)=rho2(1:nr,1,isp)-y0/nsp*rhochs(1:nr)
        enddo
C$$$C     v2 += vxc(rho2) + dv; focvxc= int (rho2*dv); dv=dvxc/drho * rhochs
C$$$      else if (lfoc .eq. 2) then
C$$$        if (ipr .gt. 40) print *, 'exchange for smooth density, foca=2:'
C$$$        if (nsp .eq. 2) then
C$$$          call dcopy(nr,rhochs,1,rhochs(1+nr),1)
C$$$        endif
C$$$c        call vxcnsp(isw,rofi,nr,rwgt,nlml,nsp,rho2,100+lxcfun,
C$$$!lm72 changed meanings of lxcfun switch.
C$$$        call vxcnsp(isw+isw2,a,rofi,nr,rwgt,nlml,nsp,rho2,10000+lxcfun,
C$$$     .  rhochs,focexc,focex,focec,focvxc,rep2,rep2x,rep2c,rmu2,v2,
C$$$     .  fl,qs)
C$$$        if (ipr .ge. 40)
C$$$     .  write(stdo,941) sumh,focexc(1)+focexc(2),focvxc(1)+focvxc(2)
C$$$  941   format(' smH core xc integrals',3f12.6)
      else
        call rxi('locpt2: cannot handle lfoc = ',lfoc)
      endif
      if (isw .eq. 1) then
        call vxcns5(1,31,'rhov*vxc',nlml,nsp,nr,rofi,rwgt,rho2,fl,xi,tmp)
        rvvxc = rvvxc - (tmp(1) + tmp(2))
        do  isp = 1, nsp
          call vxcns5(1,31,'rhov*exc',nlml,1,nr,rofi,rwgt,rho2(1,1,isp),
     .    fl(1,1,3),xi,tmp(isp))
        enddo
        rveps = rveps - (tmp(1) + tmp(2))
      endif
      if(debug) print *,'locpt2: 666661111'
      call poppr
C --- Integrals over core times effective potential ---
      vefc1 = 0d0
      do  isp = 1, nsp
        do  i = 2, nr
          ves1 = y0*v1(i,1,isp) - 2*z/rofi(i)
          vefc1 = vefc1 + rwgt(i)*rhoc(i,isp)*ves1
        enddo
      enddo
C --- Integrals involving the full nonspherical potential ---
      if(debug) print *,'locpt2: 7777'
      vefv2 = 0d0
      vefv1 = 0d0
      if (ipr .ge. 40 .and. nsp .eq. 1) write(stdo,351)
      if (ipr .ge. 40 .and. nsp .eq. 2) write(stdo,353)
      do  ilm = 1, nlml
        do  isp = 1, nsp
          rvtr(isp) = 0d0
          rvsm(isp) = 0d0
          do  i = 2, nr
            vtr = v1(i,ilm,isp)
            if (ilm .eq. 1) vtr = vtr - srfpi*2*z/rofi(i)
            rvtr(isp) = rvtr(isp) + rwgt(i)*rho1(i,ilm,isp)*vtr
            rvsm(isp) = rvsm(isp) + rwgt(i)*rho2(i,ilm,isp)*v2(i,ilm,isp)
          enddo
          vefv1 = vefv1 + rvtr(isp)
          vefv2 = vefv2 + rvsm(isp)
        enddo
        top = dmax1(dabs(rvsm(1)),dabs(rvtr(1)))
        if (ipr.ge.40 .and. top.ge.1d-6 .and. nsp .eq. 1)
     .  write(stdo,350) ilm,rvtr(1),rvsm(1)
        if (ipr.ge.40 .and. top.ge.1d-6 .and. nsp .eq. 2)
     .  write(stdo,352) ilm,rvtr(1),rvtr(2),rvtr(1)+rvtr(2),
     .  rvsm(1),rvsm(2),rvsm(1)+rvsm(2)
  350   format(i4,3x,2f15.6)
  352   format(i4,3x,3f12.6,2x,3f12.6,2x)
  351   format(/' ilm',09x,'rho*vtrue',07x,'rho*vsm')
  353   format(/' ilm',19x,'rho*vtrue',30x,'rho*vsm'/13x,
     .  'spin1',7x,'spin2',7x,'tot',11x,'spin1',7x,'spin2',7x,'tot')
      enddo
C ... Smooth xc potential includes foca head; undo in integral
      vefv2  = vefv2 - focvxc(1)-focvxc(2)
      rhoexc(1) = rep1(1) - rep2(1)
      rhoexc(2) = rep1(2) - rep2(2)
      rhoex(1)  = rep1x(1) - rep2x(1)
      rhoex(2)  = rep1x(2) - rep2x(2)
      rhoec(1)  = rep1c(1) - rep2c(1)
      rhoec(2)  = rep1c(2) - rep2c(2)
      rhovxc(1) = rmu1(1) - rmu2(1)
      rhovxc(2) = rmu1(2) - rmu2(2)
      xcore  = vefc1
      vefv2  = vefv2 + sgpotb
      valvef = vefv1 - vefv2
C --- Charges, printout ---
      if (ipr .ge. 40) then
        write(stdo,251)
        write(stdo,250) rep1(1)+rep1(2),rep2(1)+rep2(2),
     .  rhoexc(1)+rhoexc(2),
     .  rmu1(1),rmu2(1),rhovxc(1)
        if (nsp .eq. 2) write(stdo,253) rmu1(2),rmu2(2),rhovxc(2),
     .  rmu1(1)+rmu1(2),rmu2(1)+rmu2(2),
     .  rhovxc(1)+rhovxc(2)
        write(stdo,252) vefv1,vefv2,valvef,qv1,qv2,qloc
        if (nsp .eq. 2) write(stdo,254) a1,a2,aloc,alocc
        write(stdo,255) qcor1,qcor2,qlocc
  251   format(/' local terms:     true',11x,'smooth',9x,'local')
  250   format(' rhoeps:  ',3f15.6/' rhomu:   ',3f15.6)
  253   format(' spin2:   ',3f15.6/' total:   ',3f15.6)
  252   format(' val*vef  ',3f15.6/' val chg: ',3f15.6)
  254   format(' val mom: ',3f15.6,'    core:',f11.6)
  255   format(' core chg:',3f15.6)
      endif
      deallocate (fl)
      call tcx('locpt2')
      end subroutine locpt2


      subroutine elfigr(nc,stdo,z,efg1)
C- Computation of electric field gradient
C ----------------------------------------------------------------------
Ci Inputs
Ci   nc    :number of classes or sites
Ci   stdo  :standard output
Ci   z     :nuclear charge
Cio Inputs/Outputs
Cio  efg1  :input:  l=2 part of electric field at nucleus
Cio        :output: Electric field gradient
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Sep 08 Adapted from old FP (A Svane)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nc,stdo
      double precision z(nc),efg1(5,1)
C ... Local parameters
      integer i,ifi,ic,ifesn,j,n,ier,imax,ixx,iyy
      double precision v(3,3),vi(3,3),d(3),e(3),e2(3),tau(2,3)
      double precision tr(3,3),ti(3,3),da(3)
      double precision conv1,emmsfe,emmssn,Qfe,Qsn,pi,f0,s3,conv2,conv3,
     .ax,bx,cx,dx,ex,dmax,eta,split
      data conv1 /162.1d0/
      data emmsfe,emmssn /4.8d-8,7.97d-8/
      data Qfe,Qsn /0.21d-28,-0.08d-28/

      pi = 4d0*datan(1d0)
C     remember sign:
      f0 = -dsqrt(15d0/16d0/pi)
      s3 = dsqrt(3d0)
      conv2 = conv1*3d0

C      if (iprint() .lt. 30) return
      do  i = 1, 1
        ifi = stdo
        if(i .eq. 2) ifi=71
        write(ifi,'(/," Electric Field Gradient:")')
        write(ifi,'(" Site ",5x,"Tensor axes",5x,"esu/cm^2",2x,
     .  "  V/m^2  ",4x,"eta",4x,"line splitting")')
        write(ifi,'("      ",5x,"           ",5x," x10^13 ",2x,
     .  "  x10^19 ",4x,"   ",4x,"    (mm/s)")')
      enddo

      do  ic = 1, nc
c      write(ifi,*) ' ic,z=',ic,z(ic)
        if (z(ic) .gt. 0.01d0) then
          ifesn = 1
          if (dabs(z(ic)-26d0) .lt. 0.01d0) then
            conv3 = conv2*1.0d19*0.5d0*Qfe/emmsfe
          else if (dabs(z(ic)-50d0) .lt. 0.01d0) then
            conv3 = conv2*1.0d19*0.5d0*Qsn/emmssn
          else
            ifesn = 0
            conv3 = 0d0
          endif

C      ax=efg1(1,ic)+efg2(1,ic)
C      bx=efg1(2,ic)+efg2(2,ic)
C      cx=efg1(3,ic)+efg2(3,ic)
C      dx=efg1(4,ic)+efg2(4,ic)
C      ex=efg1(5,ic)+efg2(5,ic)
          ax = efg1(1,ic)
          bx = efg1(2,ic)
          cx = efg1(3,ic)
          dx = efg1(4,ic)
          ex = efg1(5,ic)
          v(1,1) = 2*ex-2/s3*cx
          v(2,2) = -2*ex-2/s3*cx
          v(3,3) = 4/s3*cx
          v(1,2) = 2*ax
          v(1,3) = 2*dx
          v(2,3) = 2*bx
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
          call htridi(n,n,v,vi,d,e,e2,tau)
          call imtql2(n,n,d,e,tr,ier)
c      write(ifi,*) ' imtql2: ier=',ier
          if (ier .gt. 0) call rx(' ELFIGR : IER ne 0')
          call htribk(n,n,v,vi,tau,n,tr,ti)
c     write(stdo,99) ic
c 99  format(/' EFG, class = ',I3)
          do  i = 1, 3
            da(i) = dabs(d(i))
          enddo
          dmax = 0d0
          imax = 0
          do  i = 1, 3
            if (da(i) .gt. dmax) then
              dmax = da(i)
              imax = i
            endif
          enddo
c  d(imax) is field gradient (Vzz)
          ixx = mod(imax,3)+1
          iyy = mod(imax+1,3)+1
          eta = 0d0
          if(dabs(d(imax)) .gt. 1.d-2) eta = dabs((d(ixx)-d(iyy))/d(imax))
c     do  i = 1, 3
c        write(stdo,98) conv1*d(i),conv2*d(i),(tr(j,i),j=1,3)
c     enddo
c     if (ifesn.eq.0) write(stdo,97) eta
          split = conv3*da(imax)*dsqrt(1d0+eta**2/3d0)
c     if (ifesn.eq.1) write(stdo,96) eta,split
c 98  format(3X,F12.4,' 10**13 esu/cm**2',3X,F12.4,' 10**19 V/m*2',
c    +         /,12X,3F12.6)
c 97  format(/' eta = ',F12.4/)
c 96  format(/' eta = ',F12.4,' line splitting = ',F12.4,' mm/s'/)

          do  i = 1, 3
            if(i .eq. 1) then
              write(ifi,'(i4,3x,3f6.3,2x,f8.2,2x,f8.2,5x,f6.3,5x,f8.5)')
     .        ic,(tr(j,i),j=1,3),conv1*d(i),conv2*d(i),eta,split
C      write(71,'(i4,3x,3f6.3,2x,f8.2,2x,f8.2,5x,f6.3,5x,f8.5)')
C     .      ic,(tr(j,i),j=1,3),conv1*d(i),conv2*d(i),eta,split
            else
              write(ifi,'(4x,3x,3f6.3,2x,f8.2,2x,f8.2,5x)')
     .        (tr(j,i),j=1,3),conv1*d(i),conv2*d(i)
C      write(71,'(4x,3x,3f6.3,2x,f8.2,2x,f8.2,5x)')
C     .         (tr(j,i),j=1,3),conv1*d(i),conv2*d(i)
            endif
          enddo
          write(ifi,'(1x)')
        endif
      enddo

      end subroutine elfigr


      end module m_locpot
