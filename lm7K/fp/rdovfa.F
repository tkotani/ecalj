Cdelw1 var ['orhofa', '1', 'orhofa(is)', '-nrmx*nsp', 'real(8)', 'rv_w_', 'rdovfa']
Cdelw1 var ['orhoca', '1', 'orhoca(is)', '-nrmx*nsp', 'real(8)', 'rv_w_', 'rdovfa']
Cdelw1 var ['ov0a', '1', 'ov0a(is)', '-nrmx*nsp', 'real(8)', 'rv_w_', 'rdovfa']
Cdelw1 var ['ov0', '1', 'ov0', 'nr*nsp', 'real(8)', 'rv_w_', 'rdovfa']
Cdelw1 var ['ov1', '1', 'ov1', 'nr*nsp', 'real(8)', 'rv_w_', 'rdovfa']
Cdelw1 var ['orwgt', '2', 'orwgt', 'nr', 'real(8)', 'rv_w_', 'rdovfa']
Cdelw1 var ['ocv', '0', 'ocv', 'ng*nsp', 'complex(8)', 'zv_w_', 'rdovfa']
Cdelw1 var ['ogv', '-1', 'undef', '0', 'undef', 'x', 'rdovfa']
Cdelw1 var ['okv', '-1', 'undef', '0', 'undef', 'x', 'rdovfa']
Cdelw1 var ['osmrho', '-1', 'undef', '0', 'undef', 'x', 'rdovfa']
Cdelw1 undel [['undef', 'wref'], 'ogv', 'rdovfa']
Cdelw1 undel [['undef', 'wref'], 'okv', 'rdovfa']
Cdelw1 undel [['undef', 'wref'], 'osmrho', 'rdovfa']
Cdelw1 undel [['call', 'array', 'wref'], 'orhoat', 'rdovfa']
Cdelw1 undel [['array', 'wref', 'call'], 'orhofa', 'rdovfa']
Cdelw1 undel [['array', 'wref', 'call'], 'orhoca', 'rdovfa']
Cdelw1 undel [['array', 'wref', 'call'], 'ov0a', 'rdovfa']
Cdelw1 undel [['wref'], 'ov0', 'rdovfa']
Cdelw1 undel [['wref'], 'ov1', 'rdovfa']
Cgetarg...info...           structure ['rdovfa', 'pot', 'spot']
Cgetarg...info...           structure ['rdovfa', 'lat', 'slat']
Cgetarg...info...           structure ['rdovfa', 'spec', 'sspec']
Cgetarg...info...           structure ['rdovfa', 'site', 'ssite']
Cgetarg...info...           use_to_add ['rdovfa', 'm_struc_def']
      subroutine rdovfa(nbas,nspec,ssite,sspec,slat,spot,qbg,orhoat)
       
       use m_struc_def  !Cgetarg

C- Read and overlap free atom densities.
C  Also allocates orhoca with free-atom core density.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   ssite :struct containing site-specific information
Ci   sspec :struct containing species-specific information
Ci   slat  :struct containing information about the lattice
Ci   qbg   :constant background charge
Co Outputs
Co   orhoat:vector of offsets containing site density, in standard
Co          3-component form (true rho, smoothed rho, core rho)
Co   smrho :smoothed interstitial density
Co         :* for smrho = smoothed mesh density, smrho is complex and
Co         :  smrho = smrho(k1,k2,k3)
Cl Local variables
Ci   k1,k2,k3 dimensions of smrho for smoothed mesh density
Cr Remarks
Cu Updates
Cu   12 May 07 package mpi-specific calls
Cu   02 Jan 06 generates core magnetic moment, checks against spec->qcor
Cu   01 Jul 05 Zero-radius empty spheres treated as local orbitals
Cu   12 Apr 03 (WRL) Added constant charge background
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   13 Jun 00 spin polarized
Cu   21 Apr 00 Adapted from nfp rdovfa.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer procid, master, mpipid
Cdelw1 do not change  [['call', 'array', 'wref'], 'orhoat', 'rdovfa']
      integer nbas,nspec,orhoat(3,nbas)
Cgetarg       double precision ssite(1),sspec(1),slat(1),spot(1),qbg
       real(8):: qbg 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat
       type(s_pot)::spot

C ... Local parameters
      integer nrmx, n0
      parameter ( nrmx=1501, n0=10 )
Cdelw1 do not change  [['array', 'wref', 'call'], 'orhofa', 'rdovfa']
Cdelw1 do not change  [['array', 'wref', 'call'], 'ov0a', 'rdovfa']
Cdelw1 do not change  [['array', 'wref', 'call'], 'orhoca', 'rdovfa']
      integer nxi(nspec),orhofa(nspec),ov0a(nspec),orhoca(nspec)
      double precision rsmfa(nspec),pnu(n0,2),exi(n0,nspec),
     .  hfc(n0,2,nspec),hfct(n0,2,nspec)
      character*8 spid(nspec),spidr
Cdelw1 do not change  [['wref'], 'ov0', 'rdovfa']
Cdelw1 do not change  [['wref'], 'ov1', 'rdovfa']
Cdelw1 do not change  [['undef', 'wref'], 'ogv', 'rdovfa']
Cdelw1 do not change  [['undef', 'wref'], 'okv', 'rdovfa']
Cdelw1 do not change  [['undef', 'wref'], 'osmrho', 'rdovfa']
Cdelw1       integer ipr,iprint,ngabc(3),n1,n2,n3,k1,k2,k3,nsp,nglob,lgunit,
Cdelw1      .  stdo,stdl,i,ifi,fopna,is,nr,lfoc,nr0,i1,nch,ib,igetss,lmxl,nlml,
Cdelw1      .  ov0,ov1,orwgt,ng,ogv,okv,ocv,iofa,osmrho,kcor,lcor
       integer:: ipr , iprint , ngabc(3) , n1 , n2 , n3 , k1 , k2 , 
     . k3 , nsp , nglob , lgunit , stdo , stdl , i , ifi , fopna , is 
     . , nr , lfoc , nr0 , i1 , nch , ib , igetss , lmxl , nlml , ov0 
     . , ov1 , ng , ogv , okv , iofa , osmrho , kcor , lcor 
       real(8) ,allocatable :: rv_w_orwgt(:)
       complex(8) ,allocatable :: zv_w_ocv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,plat(3,3),a,rmt,z,rfoc,z0,rmt0,a0,qc,ccof,
     .  ceh,stc,ztot,ctot,corm,sum,fac,sum1,sum2,sqloc,dq,vol,smom,
     .  slmom,qcor(2)
      character msg*23, strn*120
      logical mlog,cmdopt,lfail
C ... Heap
      integer w(1)
      common /w/ w

      call tcn('rdovfa')
      ipr   = iprint()
      stdo  = lgunit(1)
      stdl  = lgunit(2)
      nsp   = nglob('nsp')
      msg   = '         File mismatch:'
Cgetarg...info...           integer :: osmrho , integer(8) :: spot%osmrho 1
Cgetarg       call upack1('pot osmrho',spot,osmrho)
       
Cdelw1 warning(1) , probably  osmrho  is not defined yet at linenumber= 78
Cdelw1 do not change osmrho because of ['undef', 'wref']
       osmrho=spot%osmrho

      procid = mpipid(1)
C     numprocs = mpipid(0)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)

      if (ipr .ge. 10) write(stdo,700)
  700 format(/' rdovfa: read and overlap free-atom densities',
     .        ' (mesh density) ...')
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg       call upack('lat alat plat vol nabc',slat,alat,plat,vol,ngabc,0)
       
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       vol=slat%vol
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 

      call fftz30(n1,n2,n3,k1,k2,k3)
      call dpzero(hfc,n0*2*nspec)
      call dpzero(pnu,n0*2)

C  --- Read free-atom density for all species ---
C      if (ipr .ge. 30) write(stdo,*) ' '
       if (procid .eq. master) then
         ifi = fopna('atm',-1,0)
         rewind ifi
       endif
       do  10  is = 1, nspec
Cdelw1 do not change orhofa because of ['array', 'wref', 'call']
        call defrr(orhofa(is),   -nrmx*nsp)
Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
        call defrr(orhoca(is),   -nrmx*nsp)
Cdelw1 do not change ov0a because of ['array', 'wref', 'call']
        call defrr(ov0a(is),     -nrmx*nsp)
Cgetarg         call spacks(0,'spec name',sspec,spid(is),is,is)
         do i_spacks=is,is 
         call spacks_copy('u',sspec(i_spacks)%name,is,is,spid(is),i_spacks)
         enddo

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg         call upack('spec a nr rmt z',sspec,is,a,nr,rmt,z)
         
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt
         z=sspec(is)%z

Cgetarg...info...           integer :: lfoc , integer(8) :: sspec%lfoca 1
Cgetarg...info...           real(8) :: rfoc , real(8) :: sspec%rfoca 1
Cgetarg         call upack('spec lfoca rfoca',sspec,is,lfoc,rfoc,0,0)
         
         lfoc=sspec(is)%lfoca
         rfoc=sspec(is)%rfoca

        lfail = .false.
        if (procid .eq. master) then
        if (z .eq. 0 .and. rmt .eq. 0) then
          nxi(is) = 0
          call dpzero(exi(1,is),n0)
          call dpzero(hfc(1,1,is),2*n0)
          call dpzero(hfct(1,1,is),2*n0)
          rsmfa(is) = 0
          z0=0
          rmt0=0
          a0=0
          nr0=0
          qc=0
          ccof=0
          ceh=0
          stc=0
Cdelw1 do not change ov0a because of ['array', 'wref', 'call']
          call rlse(ov0a(is))
Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
          call rlse(orhoca(is))
Cdelw1 do not change orhofa because of ['array', 'wref', 'call']
          call rlse(orhofa(is))
        else
Cdelw1 do not change orhofa because of ['array', 'wref', 'call']
Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
Cdelw1 do not change ov0a because of ['array', 'wref', 'call']
          lfail = 
     .     (iofa(spidr,n0,nxi(is),exi(1,is),hfc(1,1,is),
     .      hfct(1,1,is),rsmfa(is),z0,rmt0,a0,nr0,qc,ccof,
     .      ceh,stc,w(orhofa(is)),w(orhoca(is)),w(ov0a(is)),ifi)
     .      .lt. 0)
        endif
        endif
        call mpibc1(nr0,1,2,mlog,'rdovfa','nr0')
        call mpibc1(lfail,1,1,mlog,'rdovfa','read error')
        if (lfail) call rxs('missing species data, species ',spid(is))
C       i = mpipid(3)
C   ... Broadcast file data
        call mpibc1(nxi(is),1,2,mlog,'rdovfa','nxi')
        call mpibc1(exi(1,is),nxi(is),4,mlog,'rdovfa','exi')
        call mpibc1(hfc(1,1,is),nsp*n0,4,mlog,'rdovfa','hfc')
        call mpibc1(hfct(1,1,is),nsp*n0,4,mlog,'rdovfa','hfct')
        call mpibc1(rsmfa(is),1,4,mlog,'rdovfa','rsmfa')
        call mpibc1(a0,1,4,mlog,'rdovfa','a0')
Cdelw1 do not change orhofa because of ['array', 'wref', 'call']
        call mpibc1(w(orhofa(is)),nr0*nsp,4,mlog,'rdovfa','rhofa')
Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
        call mpibc1(w(orhoca(is)),nr0*nsp,4,mlog,'rdovfa','rhoca')
Cdelw1 do not change ov0a because of ['array', 'wref', 'call']
        call mpibc1(w(ov0a(is)),nr0*nsp,4,mlog,'rdovfa','v0a')
        i = mpipid(3)

C ...   Defaults
        if (procid .eq. master) then
          call strip(spid(is),i1,nch)
          if (ipr .ge. 30 .and. rmt0 .ne. 0)
     .      write(stdo,400) spid(is)(1:nch),spidr,rmt0,nr0,a0
  400     format(' rdovfa: expected ',a,',',T27,' read ',a,
     .           ' with rmt=',f8.4,'  mesh',i6,f7.3)
        endif
        if (nr .le. 0)   nr = nr0
        if (a .le. 1d-6) a = a0
        if (z .eq. 0 .and. rmt .eq. 0) then
          a = 0
          nr = 0
        endif

C ...   Sanity checks
        if (procid .eq. master) then
          call fsanrg(z0,z,z,0d-9,msg,'z',.true.)
          call fsanrg(rmt0,rmt,rmt,1d-6,msg,'rmt',.true.)
          call fsanrg(a0,a,a,0d-9,msg,'a',.true.)
          call isanrg(nr0,nr,nr,msg,'nr',.true.)
        endif

C ...   Pack into sspec
Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: qc , real(8) :: sspec%qc 1
Cgetarg         call pack5('spec a nr qc',sspec,is,a,nr,qc,0)
         
         sspec(is)%a=a 
         sspec(is)%nr=nr 
         sspec(is)%qc=qc 

Cgetarg...info...           integer :: nxi(nspec) , integer(8) :: sspec%nxi 1
Cgetarg...info...           real(8) :: exi(n0,nspec) , real(8) :: sspec%exi 10
Cgetarg...info...           real(8) :: hfc(n0,2,nspec) , real(8) :: sspec%chfa 20
Cgetarg...info...           real(8) :: rsmfa(nspec) , real(8) :: sspec%rsmfa 1
Cgetarg         call pack5('spec nxi exi chfa rsmfa',sspec,is,nxi(is),
Cgetarg      .    exi(1,is),hfc(1,1,is),rsmfa(is))
         
         sspec(is)%nxi=nxi(is) 
         i_copy_size=size(sspec(is)%exi) 
         call dcopy(i_copy_size,exi(1,is),1,sspec(is)%exi,1) 
         i_copy_size=size(sspec(is)%chfa) 
         call dcopy(i_copy_size,hfc(1,1,is),1,sspec(is)%chfa,1) 
         sspec(is)%rsmfa=rsmfa(is) 

Cgetarg...info...           real(8) :: ccof , real(8) :: sspec%ctail 1
Cgetarg...info...           real(8) :: ceh , real(8) :: sspec%etail 1
Cgetarg...info...           real(8) :: stc , real(8) :: sspec%stc 1
Cgetarg...info...           integer :: orhoca(nspec) , integer(8) :: sspec%orhoc 1
Cgetarg         call pack5('spec ctail etail stc orhoc',sspec,is,ccof,ceh,stc,
Cgetarg      .    orhoca(is))
         
         sspec(is)%ctail=ccof 
         sspec(is)%etail=ceh 
         sspec(is)%stc=stc 
Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
         sspec(is)%orhoc=orhoca(is) 

   10 continue
C     Wait for all proccesses to synchronize
      i = mpipid(3)
C     Re-broadcast entire species structure, and arrays used below
      call mpibc1(sspec,nspec*nint(sspec(1)%size),4,0,'rdovfa','sspec')
      if (procid .eq. master) then
        call fclose(ifi)
      endif

C --- Define arrays for local densities rho1,rho2,rhoc and v0,v1 ---
      ztot = 0d0
      ctot = 0d0
      corm = 0d0
      do  20  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
         
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg...info...           real(8) :: qc , real(8) :: sspec%qc 1
Cgetarg         call upack('spec lmxl z qc',sspec,is,lmxl,z,qc,0)
         
         lmxl=sspec(is)%lmxl
         z=sspec(is)%z
         qc=sspec(is)%qc

Cgetarg...info...           integer :: lfoc , integer(8) :: sspec%lfoca 1
Cgetarg         call upack('spec lfoca',sspec,is,lfoc,0,0,0)
         
         lfoc=sspec(is)%lfoca

C        call upack('spec ctail etail stc orhoc',sspec,is,ccof,ceh,stc,
C     .      orca)

        nlml = (lmxl+1)**2
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
        call defrr (orhoat(1,ib),  nr*nlml*nsp)
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
        call defrr (orhoat(2,ib),  nr*nlml*nsp)
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
        call defrr (orhoat(3,ib),  nr*nsp)
Cdelw1 do not change ov0 because of ['wref']
        call defrr (ov0,           nr*nsp)
Cdelw1 do not change ov1 because of ['wref']
        call defrr (ov1,           nr*nsp)

C       Core magnetic moment (possible if magnetized core hole)
        if (nsp .eq. 2 .and. lmxl .gt. -1) then
Cdelw1           call defrr(orwgt, nr)
           allocate(rv_w_orwgt(nr))
           if (nr<0) rv_w_orwgt(:)=0.0d0

Cdelw1           call radwgt(rmt,a,nr,w(orwgt))
           call radwgt ( rmt , a , nr , rv_w_orwgt ) 

Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
Cdelw1           call radsum(nr,nr,1,nsp,w(orwgt),w(orhoca(is)),sum)
           call radsum ( nr , nr , 1 , nsp , rv_w_orwgt , w ( orhoca ( is 
     .     ) ) , sum ) 

Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
Cdelw1           call radsum(nr,nr,1,1,w(orwgt),w(orhoca(is)),sum1)
           call radsum ( nr , nr , 1 , 1 , rv_w_orwgt , w ( orhoca ( is 
     .     ) ) , sum1 ) 

          sum2 = sum - sum1
          call gtpcor(sspec,is,kcor,lcor,qcor)

          if (dabs(qcor(2)-(sum1-sum2)) .gt. 0.01d0) then
            call info5(10,0,0,' (warning) core moment mismatch spec %i:'
     .        //'  input file=%;6d  atom file=%;6d',
     .        is,qcor(2),sum1-sum2,0,0)
          endif
          corm = corm + qcor(2)
Cdelw1 rlse name= orwgt old_list= orhofa orhoca ov0a orhoat ov0 ov1 orwgt 
Cdelw1 rlse name= orwgt new_list= orhofa orhoca ov0a orhoat ov0 ov1 
Cdelw1           call rlse(orwgt)
           if (allocated(rv_w_orwgt)) deallocate(rv_w_orwgt)

        endif

        if (lmxl .gt. -1) then
Cdelw1 do not change ov0a because of ['array', 'wref', 'call']
Cdelw1 do not change ov0 because of ['wref']
          call dpcopy(w(ov0a(is)),  w(ov0),  1,nr*nsp,1d0)
Cdelw1 do not change ov0a because of ['array', 'wref', 'call']
Cdelw1 do not change ov1 because of ['wref']
          call dpcopy(w(ov0a(is)),  w(ov1),  1,nr*nsp,1d0)
Cdelw1 3 is digit, dropped
Cdelw1 do not change orhoca because of ['array', 'wref', 'call']
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
          call dpcopy(w(orhoca(is)),w(orhoat(3,ib)),1,nr*nsp,1d0)
  
          if (lfoc .eq. 0) then
Cdelw1             call defrr(orwgt, nr)
             allocate(rv_w_orwgt(nr))
             if (nr<0) rv_w_orwgt(:)=0.0d0

Cdelw1             call radwgt(rmt,a,nr,w(orwgt))
             call radwgt ( rmt , a , nr , rv_w_orwgt ) 

Cdelw1 3 is digit, dropped
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
Cdelw1             call radsum(nr,nr,1,nsp,w(orwgt),w(orhoat(3,ib)),sum)
             call radsum ( nr , nr , 1 , nsp , rv_w_orwgt , w ( orhoat ( 3 
     .       , ib ) ) , sum ) 

            fac = 1d0
            if(dabs(sum) .gt. 1d-7) fac = qc/sum
            if (ipr .ge. 40) write(stdo,787) is,qc,sum,fac
  787       format(' scale foca=0 core species',i2,': qc,sum,scale=',
     .        3f12.6,f12.6)
Cdelw1 3 is digit, dropped
Cdelw1 3 is digit, dropped
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
            call dpcopy(w(orhoat(3,ib)),w(orhoat(3,ib)),1,nr*nsp,fac)
Cdelw1 rlse name= orwgt old_list= orhofa orhoca ov0a orhoat ov0 ov1 orwgt 
Cdelw1 rlse name= orwgt new_list= orhofa orhoca ov0a orhoat ov0 ov1 
Cdelw1             call rlse(orwgt)
             if (allocated(rv_w_orwgt)) deallocate(rv_w_orwgt)

          endif
        endif

Cgetarg...info...           integer :: ov0 , integer(8) :: ssite%ov0 1
Cgetarg...info...           integer :: ov1 , integer(8) :: ssite%ov1 1
Cgetarg         call pack5('site ov0 ov1',ssite,ib,ov0,ov1,0,0)
         
Cdelw1 do not change ov0 because of ['wref']
         ssite(ib)%ov0=ov0 
Cdelw1 do not change ov1 because of ['wref']
         ssite(ib)%ov1=ov1 

        ztot = ztot+z
        ctot = ctot+qc

C     end loop over sites
   20 continue

C --- Overlap smooth hankels to get smooth interstitial density ---
Cgetarg...info...           integer :: ng , integer(8) :: slat%ng 1
Cgetarg...info...           integer :: ogv , integer(8) :: slat%ogv 1
Cgetarg...info...           integer :: okv , integer(8) :: slat%okv 1
Cgetarg       call upack('lat ng ogv okv',slat,ng,ogv,okv,0,0)
       
       ng=slat%ng
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 347
Cdelw1 do not change ogv because of ['undef', 'wref']
       ogv=slat%ogv
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 348
Cdelw1 do not change okv because of ['undef', 'wref']
       okv=slat%okv

Cdelw1       call defcc(ocv, ng*nsp)
       allocate(zv_w_ocv(ng*nsp))
       if (ng*nsp<0) zv_w_ocv(:)=0.0d0

Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 351
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 351
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1       call ovlpfa(ssite,slat,nbas,nxi,n0,exi,hfc,rsmfa,ng,ng,w(ogv),
Cdelw1      .  w(ocv))
       call ovlpfa ( ssite , slat , nbas , nxi , n0 , exi , hfc , rsmfa 
     . , ng , ng , w ( ogv ) , zv_w_ocv ) 


ccccccccccccccc
c        is=1
c        call upack('spec nr',sspec,is,nr,0,0,0)
c        print *,' qqq: aaaaaaaaaaaaaaaaaaaaaaa is nr='
c     &     ,is, nr,ng,nsp,k1,k2,k3,okv,ocv,osmrho
c        call wkchk('qqq: aaa')
c        print *,' qqq: okv ocv osmrho=', okv,ocv,osmrho
ccccccccccccccc

Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 363
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  osmrho  is not defined yet at linenumber= 363
Cdelw1 do not change osmrho because of ['undef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 363
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  osmrho  is not defined yet at linenumber= 363
Cdelw1 do not change osmrho because of ['undef', 'wref']
Cdelw1       call gvputf(ng,nsp,w(okv),k1,k2,k3,w(ocv),w(osmrho))
       call gvputf ( ng , nsp , w ( okv ) , k1 , k2 , k3 , zv_w_ocv 
     . , w ( osmrho ) ) 


ccccccccccccccc
c        is=1
c        call upack('spec nr',sspec,is,nr,0,0,0)
c        call wkchk('qqq: 1bb')
c        print *,'qqq: 2 aaaaaaaaaaaaaaaaaaaaaaa is nr=',is, nr
c        stop 'xxxxxxx'
ccccccccccccccc
Cdelw1 rlse name= ocv old_list= orhofa orhoca ov0a orhoat ov0 ov1 ocv 
Cdelw1 rlse name= ocv new_list= orhofa orhoca ov0a orhoat ov0 ov1 
Cdelw1       call rlse(ocv)
       if (allocated(zv_w_ocv)) deallocate(zv_w_ocv)




C ... FFT to real-space mesh
Cdelw1 warning(1) , probably  osmrho  is not defined yet at linenumber= 377
Cdelw1 do not change osmrho because of ['undef', 'wref']
      call fftz3(w(osmrho),n1,n2,n3,k1,k2,k3,nsp,0,1)
C ... Add compensating uniform electron density to compensate background
Cdelw1 warning(1) , probably  osmrho  is not defined yet at linenumber= 379
Cdelw1 do not change osmrho because of ['undef', 'wref']
      call addbkgsm(w(osmrho),k1,k2,k3,nsp,qbg,vol,-1d0)
C ... integrate
Cdelw1 warning(1) , probably  osmrho  is not defined yet at linenumber= 381
Cdelw1 do not change osmrho because of ['undef', 'wref']
      call mshint(vol,nsp,n1,n2,n3,k1,k2,k3,w(osmrho),sum1,sum2)
      if (nsp .eq. 2) then
Cdelw1 warning(1) , probably  osmrho  is not defined yet at linenumber= 383
Cdelw1 do not change osmrho because of ['undef', 'wref']
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,w(osmrho),smom,sum2)
        smom = 2*smom - sum1
      endif
C --- Set up local densities using rmt from atm file ---


Cdelw1 do not change orhofa because of ['array', 'wref', 'call']
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
      call ovlocr(nbas,ssite,sspec,slat,n0,nxi,exi,hfc,rsmfa,orhofa,
     .  orhoat,sqloc,slmom)
C --- Add compensating uniform electron density to compensate background
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
      call adbkql(orhoat,nbas,nsp,qbg,vol,-1d0,sspec,ssite)
      if (abs(qbg) .ne. 0) call info(10,0,0,' Uniform '//
     .  'density added to neutralize background, q=%;6,6d',qbg,0)

C --- Print charges ---
      dq = sum1+sqloc+ctot-ztot+qbg
      if (nsp .eq. 1) then
        if (ipr .ge. 10)
     .  write(stdo,895) sum1,sqloc,sum1+sqloc,ctot,-ztot,qbg,dq
  895   format(/' Smooth charge on mesh:    ',f16.6
     .     /    ' Sum of local charges:     ',f16.6
     .     /    ' Total valence charge:     ',f16.6
     .     /    ' Sum of core charges:      ',f16.6
     .     /    ' Sum of nuclear charges:   ',f16.6
     .     /    ' Homogeneous background:   ',f16.6
     .     /    ' Deviation from neutrality:',f16.6)
        if (ipr .ge. 10) write (stdl,710) sum1+sqloc,sum1,sqloc,qbg,dq
  710   format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .    '   bg',f10.6,'  dQ',f10.6)
      else
        if (ipr .ge. 10) write(stdo,896) sum1,smom,sqloc,slmom,
     .    sum1+sqloc,smom+slmom,ctot,corm,-ztot,qbg,dq
  896   format(/' Smooth charge on mesh:    ',f16.6,4x,'moment', f12.6,
     .     /    ' Sum of local charges:     ',f16.6,4x,'moments',f11.6,
     .     /    ' Total valence charge:     ',f16.6,4x,'moment', f12.6,
     .     /    ' Sum of core charges:      ',f16.6,4x,'moment', f12.6,
     .     /    ' Sum of nuclear charges:   ',f16.6
     .     /    ' Homogeneous background:   ',f16.6
     .     /    ' Deviation from neutrality:',f16.6)
        if (ipr .ge. 10)
     .  write (stdl,711) sum1+sqloc,sum1,sqloc,qbg,smom+slmom
  711   format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .    '   bg',f11.6,' mm',f11.6)
      endif

      if (dabs(dq) .gt. 1d-4 .and. ipr .gt. 0)
     .  call awrit1(' rdovfa (warning) overlapped'
     .  //' density not neutral'//', dq= %d',' ',80,stdo,dq)

      call tcx('rdovfa')

Cgetarg       end
Cdelw1 w_varlist remains: orhofa orhoca ov0a orhoat ov0 ov1 
Cdelw1 w_varlistundel: [ogv okv osmrho orhoat orhofa orhoca ov0a ov0 ov1]
Cdelw1 w_varlist (undel), remains: [orhofa orhoca ov0a orhoat ov0 ov1]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate ov1 because of [wref]
Cdelw1 not deallocate ov0 because of [wref]
Cdelw1 not deallocate orhoat because of [call array wref]
Cdelw1 not deallocate ov0a because of [array wref call]
Cdelw1 not deallocate orhoca because of [array wref call]
Cdelw1 not deallocate orhofa because of [array wref call]
       end subroutine rdovfa 


