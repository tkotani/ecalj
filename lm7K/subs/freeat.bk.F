Cdelw1 var ['og', '2', 'og', 'nr*2', 'real(8)', 'rv_w_', 'freats']
Cdelw1 var ['opsi', '2', 'opsi', 'nr*(lmxa+1)*nsp', 'real(8)', 'rv_w_', 'freats']
Cdelw1 var ['oh', '1', 'oh', 'nr', 'real(8)', 'rv_w_', 'optfab']
Cdelw1 var ['og', '1', 'og', '2*nr', 'real(8)', 'rv_w_', 'optfab']
Cdelw1 var ['ogp', '1', 'ogp', '2*nr*4', 'real(8)', 'rv_w_', 'optfab']
Cdelw1 var ['opsi', '1', 'opsi', 'nr*lmxa', 'real(8)', 'rv_w_', 'optfab']
Cdelw1 var ['oh', '1', 'oh', 'nr', 'real(8)', 'rv_w_', 'ftfalo']
Cdelw1 var ['og', '1', 'og', '2*nr', 'real(8)', 'rv_w_', 'ftfalo']
Cdelw1 var ['ogp', '1', 'ogp', '2*nr*4', 'real(8)', 'rv_w_', 'ftfalo']
Cdelw1 var ['opsi', '1', 'opsi', 'nr*lmxa', 'real(8)', 'rv_w_', 'ftfalo']
Cgetarg...info...           structure ['freeat', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['freeat', 'spec', 'sspec']
Cgetarg...info...           structure ['freeat', 'ham', 'sham']
Cgetarg...info...           use_to_add ['freeat', 'm_struc_def']
       subroutine freeat(sctrl,sspec,sham,spot)
        
        use m_struc_def  !Cgetarg

C- For each species, makes free atom self-consistent
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: lxcf smalit
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: rsmfa rfoca coreq z rmt a nr p q idmod lmxa pz eref 
Ci                rs3 eh3 vmtz
Ci     Stored:    name coreh z rmt a nr norp ntorb orbp
Ci     Passed to: spacks gtpcor ioorbp
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: seref
Ci   spot  :struct containing information about the potential
Co Outputs
Cl Local variables
Cl   ccof  :coefficient to fit of core tail to smoothed Hankel
Cl   ceh   :energy of core tail to smoothed Hankel
Cl   sumtc :core kinetic energy
Cr Remarks
Cu Updates
Cu   01 Feb 06 Enables renormalized free atom density
Cu   01 Jul 05 Skips spheres with Z=0 and R=0
Cu   21 Jun 04 Added fit of sm. Hankel tails to local orbitals
Cu   18 Sep 03 (ATP) Enabled partial core occupation
Cu   06 Sep 03 Constrain rsm in fit to FA wave function
Cu   18 Mar 03 Altered sign of magnetic moment to conform to std
Cu   19 Apr 02 Redesigned freats call to avoid the use of structures
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   22 Mar 01 Added printout of reference energy
Cu   10 Jun 00 spin polarized
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
Cgetarg       double precision sctrl(1),sham(1),spot(1),sspec(1)
       type(s_ctrl)::sctrl
       type(s_ham)::sham
       type(s_pot)::spot
       type(s_spec)::sspec(*)

C ... Local parameters
      integer fopna,ifi,iprint,is,lgunit,nglob,nr,nrmt,nrmx,nsp,nspec,
     .  n0,nkap0,nxi,nxi0,stdo,nrmix(2),lxcfun,igets,lmxa,kcor,lcor
      character*8 spid,chole*8
      parameter ( nrmx=1501, nxi0=10, n0=10, nkap0=3)
      double precision qc,ccof,ceh,z,rmt,rfoca,rsmfa,qcor(2),a,sumec,
     .  sumtc,eref,seref,dgets,dgetss,etot
      double precision hfc(nxi0,2),exi(nxi0),hfct(nxi0,2)
      double precision v(nrmx*2),rho(nrmx*2),rhoc(nrmx*2),rofi(nrmx*2)
      double precision pnu(n0,2),pz(n0,2),qat(n0,2),rcfa(2)
      double precision rtab(n0,2),etab(n0,2),orbp(n0,2,nkap0)
      double precision rs3,eh3,vmtz
      integer idmod(n0)
ctakao
      logical ixx
      ifi = fopna('atm',-1,0)
      rewind ifi
      nspec = nglob('nspec')
      nsp   = nglob('nsp')
      stdo = lgunit(1)
      exi(1) = -1
      exi(2) = -2
      exi(3) = -4
      exi(4) = -6
      exi(5) = -9
      exi(6) = -15
      nxi = 6
      call dpzero(hfct,2*nxi0)

      do  is = 1, nspec

C       Unpack sctrl data
Cgetarg...info...           integer :: nrmix(2) , integer(8) :: sctrl%smalit 2
Cgetarg         call upack('ctrl smalit',sctrl,nrmix,0,0,0,0)
         
         i_copy_size=size(sctrl%smalit) 
         call i8icopy(i_copy_size,sctrl%smalit,1,nrmix,1) 

Cgetarg         lxcfun = igets('ctrl lxcf',sctrl)
         lxcfun = int(sctrl%lxcf) 


C       Unpack species data
Cgetarg         call spacks(0,'spec name',sspec,spid,is,is)
         do i_spacks=is,is 
         call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
         enddo

Cgetarg...info...           real(8) :: rsmfa , real(8) :: sspec%rsmfa 1
Cgetarg...info...           real(8) :: rfoca , real(8) :: sspec%rfoca 1
Cgetarg...info...           real(8) :: qcor(2) , real(8) :: sspec%coreq 2
Cgetarg         call upack('spec rsmfa rfoca coreq',sspec,is,rsmfa,rfoca,qcor,0)
         
         rsmfa=sspec(is)%rsmfa
         rfoca=sspec(is)%rfoca
         i_copy_size=size(sspec(is)%coreq) 
         call dcopy(i_copy_size,sspec(is)%coreq,1,qcor,1) 

Cgetarg         call spacks(0,'spec coreh',sspec,chole,is,is)
         do i_spacks=is,is 
         call spacks_copy('u',sspec(i_spacks)%coreh,is,is,chole,i_spacks)
         enddo

        call gtpcor(sspec,is,kcor,lcor,qcor)
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nrmt , integer(8) :: sspec%nr 1
Cgetarg         call upack('spec z rmt a nr',sspec,is,z,rmt,a,nrmt)
         
         z=sspec(is)%z
         rmt=sspec(is)%rmt
         a=sspec(is)%a
         nrmt=sspec(is)%nr

        if (z .eq. 0 .and. rmt .eq. 0) goto 10
        call dpzero(qat,n0*2)
Cgetarg...info...           real(8) :: pnu(n0,2) , real(8) :: sspec%p 20
Cgetarg...info...           real(8) :: qat(n0,2) , real(8) :: sspec%q 20
Cgetarg...info...           integer :: idmod(n0) , integer(8) :: sspec%idmod 10
Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg         call upack('spec p q idmod lmxa',sspec,is,pnu,qat,idmod,lmxa)
         
         i_copy_size=size(sspec(is)%p) 
         call dcopy(i_copy_size,sspec(is)%p,1,pnu,1) 
         i_copy_size=size(sspec(is)%q) 
         call dcopy(i_copy_size,sspec(is)%q,1,qat,1) 
         i_copy_size=size(sspec(is)%idmod) 
         call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1) 
         lmxa=sspec(is)%lmxa

Cgetarg...info...           real(8) :: pz(n0,2) , real(8) :: sspec%pz 20
Cgetarg         call upack2('spec pz',sspec,is,pz)
         
         i_copy_size=size(sspec(is)%pz) 
         call dcopy(i_copy_size,sspec(is)%pz,1,pz,1) 

Cgetarg         eref = dgetss('spec eref',is,sspec)
         eref = (sspec(is)%eref) 

Cgetarg...info...           real(8) :: rs3 , real(8) :: sspec%rs3 1
Cgetarg...info...           real(8) :: eh3 , real(8) :: sspec%eh3 1
Cgetarg...info...           real(8) :: vmtz , real(8) :: sspec%vmtz 1
Cgetarg...info...           real(8) :: rcfa(2) , real(8) :: sspec%rcfa 2
Cgetarg         call upack('spec rs3 eh3 vmtz rcfa',sspec,is,rs3,eh3,vmtz,rcfa)
         
         rs3=sspec(is)%rs3
         eh3=sspec(is)%eh3
         vmtz=sspec(is)%vmtz
         i_copy_size=size(sspec(is)%rcfa) 
         call dcopy(i_copy_size,sspec(is)%rcfa,1,rcfa,1) 


        call freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,kcor,lcor,qcor,
     .    nrmix,1,lxcfun,z,rmt,a,nrmt,pnu,pz,qat,rs3,eh3,vmtz,rcfa,
     .    idmod,lmxa,eref,rtab,etab,hfc,hfct,nr,rofi,rho,rhoc,qc,ccof,
     .    ceh,sumec,sumtc,v,etot)

C       Repack mesh parms, in case a mesh was selected by freats
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nrmt , integer(8) :: sspec%nr 1
Cgetarg         call pack5('spec z rmt a nr',sspec,is,z,rmt,a,nrmt)
         
         sspec(is)%z=z 
         sspec(is)%rmt=rmt 
         sspec(is)%a=a 
         sspec(is)%nr=nrmt 


C   ... Pack the basis into sspec
        call dpzero(orbp,n0*2*nkap0)
C       call upack2('spec orbp lmxb',sspec,is,orbp)
        call dcopy(n0,rtab,1,orbp(1,1,1),1)
        call dcopy(n0,etab,1,orbp(1,2,1),1)
Cgetarg...info...           integer :: n0 , integer(8) :: sspec%ntorb 1
Cgetarg         call pack5('spec norp ntorb',sspec,is,2,n0,0,0)
         
         sspec(is)%norp=2 
         sspec(is)%ntorb=n0 

Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: sspec%orbp 60
Cgetarg         call pack5('spec orbp',sspec,is,orbp,0,0,0)
         
         i_copy_size=size(sspec(is)%orbp) 
         call dcopy(i_copy_size,orbp,1,sspec(is)%orbp,1) 


C   --- File write ---
        if (iprint() .gt. 40) write(stdo,230) spid
  230   format(/' write free atom data for species  ',a)
C       Copy second spin channel of rho,v
        if (nsp .eq. 2 .and. nr .gt. nrmt) then
          call dcopy(nrmt,rho(1+nr),1,rho(1+nrmt),1)
          call dcopy(nrmt,rhoc(1+nr),1,rhoc(1+nrmt),1)
          call dcopy(nrmt,v(1+nr),1,v(1+nrmt),1)
        endif
        call iofa(spid,nxi0,nxi,exi,hfc,hfct,rsmfa,z,rmt,
     .     a,nrmt,qc,ccof,ceh,sumtc,rho,rhoc,v,-ifi)

   10 continue
      enddo

cccccccccccccccccccccc
ctakao
      call fclose(ifi)
c      ixx= ioorbp(1,2,1,nspec,sspec,0,-ifi)
c      call fclose(ifi)
ccccccccccccccccccccc
c takao
      ifi = fopna('mtopara',-1,0)
      call writebasis(nspec,sspec,ifi)
      call fclose(ifi)
cccccccccccccc
      if (iprint() .gt. 30)  then
Cgetarg         seref = dgets('ham seref',sham)
         seref = sham%seref 

        call awrit1('%x%N Sum of reference energies: %1;6d',' ',80,stdo,
     .    seref)
      endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine freeat 


      subroutine freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,kcor,lcor,qcor,
     .  nrmix,lwf,lxcfun,z,rmt,a,nrmt,pnu,pz,qat,rs3,eh3,vmtz,rcfa,
     .  idmod,lmxa,eref,rtab,etab,hfc,hfct,nr,rofi,rho,rhoc,qc,ccof,ceh,
     .  sec,stc,v,etot)

C- Makes one free atom self-consistent, fits rho tails to smoothed Hankels
C ----------------------------------------------------------------------
Ci Inputs
Ci   spid  :species label (for printout)
Ci   is    :species index
Ci   nxi0  :nxi0: leading dimension of hfc,hfct
Ci   nxi   :number of hankel functions used in fitting of tails
Ci   exi   :hankel energies used in fitting of tails
Ci   rfoca :smoothing radius for hankel fit to core
Ci   rsmfa :smoothing radius for hankel fit to valence
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   nrmix :nrmix(1) = maximum number of interations in sphere
Ci         :           before giving up on self-consistency
Ci         :nrmix(2) = no prior iterations Anderson mixing in
Ci         :           self-consistency cycle.
Ci   lwf   :1 print information about wave functions
Ci   lxcfun:selects local exchange-correlation functional
Ci   z     :nuclear charge
Ci   rmt   :augmentation radius, in a.u.
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nrmt  :number of mesh points from origin to rmt
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pz    :boundary conditions for local orbitals
Ci   qat   :valence charges for each l channel
Ci   rs3   :minimum allowed smoothing radius in attaching Hankel tails
Ci         :to local orbitals
Ci   eh3   :Hankel energy when attaching Hankel tails to high-lying
Ci         :local orbitals
Ci   vmtz  :parameter used in attaching Hankel tails to local orbitals
Ci         :It is used as a constant shift to Hankel energies for the 
Ci         :fitting of local orbitals to Hankel tails. Thus vmtz
Ci         :is an estimate for the potential at the MT radius.
Ci   idmod :0,1 or 2, specifing how the enu is set for an l-channel
Ci   lmxa  :augmentation l-cutoff
Ci   eref  :reference energy (used for printout)
Co Outputs
Co   rtab  :smoothing radius for optimized wave function
Co   etab  :energy for optimized wave function
Co   hfc   :fit coeffs for valence density,
Co   hfct  :contains fit coeffs for full density (not calc. now)
Co   nr    :number of radial mesh points for spherical rho
Co   rofi  :rofi(1..nr)=radial mesh for points
Co         :rofi(nr+1..2*nr)=radial mesh weights
Co   rho   :free-atom valence density
Co   rhoc  :free-atom core density
Co   qc    :Sphere core charge
Co   ccof  :coefficient to fit of core tail to smoothed Hankel
Co   ceh   :energy of core tail to smoothed Hankel
Co   sec   :sum of core eigenvalues
Co   stc   :core kinetic energy
Co   v     :spherical potential
Cl Local variables
Cl   itab  :itab(l+1)=1  a wave function was optimzed for this l
Cl   pnul  :EITHER : pnu for valence state, OR
Cl         :local orbital if DEEPER than valence (pz<pnu)
Cr Remarks
Cu Updates
Cu   01 Feb 06 Enables renormalized free atom density
Cu   19 Apr 02 Redesigned input to avoid the use of structures
Cu   10 Apr 02 Redimensionsed etab,rtab to accomodate larger lmax
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nrmx,nrmt,is,nxi0,nxi,nrmix(2),lwf,lxcfun,n0,kcor,lcor
      parameter (nrmx=1501,n0=10)
      character*8 spid
      double precision rsmfa,rfoca,qc,ccof,ceh,sec,stc,z,rmt,a,eref,
     .  v(nrmx*2),rho(nrmx*2),rhoc(nrmx*2),hfc(nxi0,1),hfct(nxi0,1),
     .  exi(1),rtab(n0,2),etab(n0,2),rofi(nrmx*2),rs3,eh3,vmtz,qcor(2)
C ... Local parameters
      logical cmdopt
      integer ncmx,nvmx
      parameter (ncmx=50, nvmx=20)
      integer idmod(n0)
C     integer idmoz(n0)
      character str*8,strn*32
      double precision rmax,b,etot,dq,
     .  ec(ncmx),ev(nvmx),sumev,vrmax(2),exrmax(2),ekin,utot,rhoeps,
     .  amgm,rhrmx,qvt,qtot,qct,qvin,qcin,r,wt0,wt1,qtt,qtin,pnul,
     .  pzl,pnu(n0,2),qat(n0,2),pl(n0,2),rhoin(nrmx*2),
     .  rhot(nrmx*2),ql(3,n0,2),pz(n0,2),rcfa(2)
C     double precision qz(n0,2)
      integer fopna,i,ifi,ipr,iprint,isp,isw,l,lfrz,
     .  lgrad,lgunit,lmxa,lplfa,lrel,nglob,nitmax,nmix,nr,
     .  nsp,stdl,stdo,lplawv,irchan(n0)
Cdelw1       integer og,opsi
       real(8) ,allocatable :: rv_w_og(:)
       real(8) ,allocatable :: rv_w_opsi(:)

C     Parameters for optimized free-atom wave functions
      integer itab(n0,2)
C Heap
      integer w(1)
      common /w/ w
ctakao
      integer:: ipl, ipz,iplx
      real(8):: qcc

      stdo   = lgunit(1)
      stdl   = lgunit(2)
      ipr    = iprint()
      lrel   = nglob('lrel')
      nsp    = nglob('nsp')
      lfrz   = 0
      lgrad  = lxcfun/10

ctakao
      if(ipr>29) then
      print *
      do i = 1, nsp
      do l = 0, lmxa
      write(6,"(a,2i3,2f10.3)")'ttt: pnu qat=',i,l,pnu(l+1,i),qat(l+1,i)
      enddo
      enddo
      endif

C --- Get species data ----
      call dpzero(pl,2*n0)
      call dpzero(ql,2*3*n0)
      do  10  i = 1, nsp
      do  10  l = 0, lmxa
        pnul = pnu(l+1,i)
        pzl  = mod(pz(l+1,1),10d0)
        if (pzl .ne. 0) then
C         Make pz the valence state: hardwire qat
ctakao
          if( int(pzl)<=l .and. pzl/=0d0) then
             print *,' l  pz=',l,pzl
             call rx(' Int(PZ) > angular momentum l should be.')
          elseif (int(pnul-1) .eq. int(pzl)) then
ctakao
c          if (int(pnul-1) .eq. int(pzl)) then
c            if (qat(l+1,i) .ne. 0) call fexit3(-1,111,
c     .        ' Exit -1 freeat, l=%i:  nonzero charge Q=%d not allowed'
c     .        //'  with sc PZ=%;3d',l,qat(l+1,i),pzl)
ctakao
c            pnul = pzl
c            qat(l+1,1) = 4*l+2
c            if(i==1) qat(l+1,1) = 4*l+2 + qat(l+1,1)
c            qat(l+1,2) = 0
          elseif (int(pnul+1) .ne. int(pzl)) then
            call fexit3(-1,111,' Exit -1 freeat, l=%i:  '//
     .        'sc PZ=%d incompatible with valence P=%;3d',l,pzl,pnul)
          endif

        endif
        pl(l+1,i) = int(pnul) + .5d0
        ql(1,l+1,i) = qat(l+1,1)/nsp
        if (nsp .eq. 2) then
          if (pnul .eq. 0) pl(l+1,i) = pl(l+1,1)
          ql(1,l+1,i) = qat(l+1,1)/nsp - qat(l+1,2)/2*dble(2*i-3)
        endif
        ql(2,l+1,i) = 0d0
        ql(3,l+1,i) = 0d0
   10 continue

ctakao
      print *
      write(6,"(a,a)") 'SPEC_ATOM= '//trim(spid),
     & ' --- Table for atomic configuration --- '
      write(6,"(a)")' isp  l   int(P) int(PZ)   QVAL-QZ  QCORE+QZ   QZ'
     &              //'       QVAL    QCORE '
      qtot=0d0
      qc  = 0d0
      amgm= 0d0
      qzz=0d0
      do i = 1, nsp
      do l = 0, lmxa
        ipl = int(pl(l+1,i))
        ipz = int(mod(pz(l+1,1),10d0))
        qcc= (ipl-l-1)*(4*l+2)/nsp
        qc = qc+qcc
        qzz=0d0
        if(ipz/=0.and.ipz<ipl) qzz=(4*l+2)/nsp
        write(6,"(i4,i3,3x,i5,i7,f11.3,4f9.3)") i,l,ipl,ipz
     &   ,ql(1,l+1,i), qcc, qzz
     &   ,ql(1,l+1,i)+QZ, qcc-qzz
        qtot=qtot +         qcc+ql(1,l+1,i)
        amgm=amgm +(3-2*i)*(qcc+ql(1,l+1,i))
      enddo
        print * 
      enddo
c      call getqvc(nsp,n0,lmxa,z,pl,ql,0,0,kcor,lcor,qcor,
c     .            qc,qtot,amgm,w,w)
      if (ipr .ge. 20) then
        call awrit6('%N Species '//spid//'%a:  Z=%d'//
     .    '  Qc=%d  R=%1,6;6d  Q=%1;6d%?#n==2#  mom=%1;5d#%0d#%a',
     .    ' ',80,stdo,z,qc,rmt,qtot,nsp,amgm)
      endif


C --- Set up radial mesh for free atom ---
      rmax = 50d0
      if (z .lt. 10) rmax = 25
      if (z .le.  6) rmax = 20
      call pshpr(0)
      call rmesh(z,rmt,lrel,lgrad,nrmx,a,nrmt)
      call poppr
      b = rmt/(dexp(a*nrmt-a)-1d0)
      nr = 1d0+dlog(1d0+rmax/b)/a
      if (mod(nr,2).eq.0) nr = nr-1
      rmax = b*(dexp(a*(nr-1))-1d0)
      call info5(21,0,0,' mesh:   rmt=%,6;6d  rmax=%,6;6d'//
     .  '  a=%d  nr=%i  nr(rmax)=%i',rmt,rmax,a,nrmt,nr)

C --- Make atom self-consistent ---
      nitmax = nrmix(1)
      nmix = nrmix(2)
      nmix = -30
C     call pshpr(min(iprint(),40))
      ec(1) = 0
      call atomsc(.false.,n0,nsp,lmxa,z,0d0,kcor,lcor,qcor,rmax,a,nr,
     .  rofi,ec,ev,pl,ql,idmod,v,0d0,rhoin,rho,rhoc,nmix,qc,sec,stc,
     .  sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,dq,exrmax,'gue',
     .  nitmax,lfrz)
C     call poppr

      if (ipr .ge. 20)
     .  call awrit6('%?#n>=30#%N## sumev=%,6;6d  etot=%,6;6d'//
     .  '  eref=%,6;6d%?#n#  diff= %,6;6d',' ',80,stdo,ipr,sumev,etot,
     .  eref,isw(eref.ne.0),etot-eref)

      call dcopy(lmxa+1,ql,3,qat,1)
      call dcopy(lmxa+1,ql(1,1,nsp),3,qat(1,2),1)
      call awrit4('fa  Pl %n:-1d  Ql %n:-1d',
     .  ' ',80,stdl,lmxa+1,pl,lmxa+1,qat)
      if (nsp .eq. 2)
     .call awrit4('fa  Pl2%n:-1d  Ql2%n:-1d',
     .  ' ',80,stdl,lmxa+1,pl(1,nsp),lmxa+1,qat(1,nsp))

      if (dabs(dq) .gt. 1d-5 .and. iprint() .ge. 10)
     .  call awrit1(' freeat (warning) atom not neutral, Q=%d',
     .  ' ',80,lgunit(1),dq)
C .. Subtract core from density to make valence density
      do  isp = 1, nsp
      do  i = 1, nr
        rhot(i+(isp-1)*nr) = rho(i+(isp-1)*nr)
        rho(i+(isp-1)*nr) = rho(i+(isp-1)*nr)-rhoc(i+(isp-1)*nr)
      enddo
      enddo

C --- Renormalize atom density or potential ---
      call ivset(irchan,1,n0,0)
      call rnatm(pl,qat,n0,irchan,lmxa,z,a,b,rofi,ev,nr,rcfa,nsp,v,rho)
C     call prrmsh('starting total rho',rofi,rhot,nr,nr,nsp)
      do  isp = 1, nsp
      do  i = 1, nr
        rhot(i+(isp-1)*nr) = rho(i+(isp-1)*nr)+rhoc(i+(isp-1)*nr)
      enddo
      enddo
C     call prrmsh('ending total rho',rofi,rhot,nr,nr,nsp)

C --- Print info about free-atom wavefunctions ---
      if (lwf .ne. 0) then
      if (ipr .gt. 30) then
Cdelw1       call defrr(og,  nr*2)
       allocate(rv_w_og(nr*2))
       if (nr*2<0) rv_w_og(:)=0.0d0

Cdelw1       call defrr(opsi, nr*(lmxa+1)*nsp)
       allocate(rv_w_opsi(nr*(lmxa+1)*nsp))
       if (nr*(lmxa+1)*nsp<0) rv_w_opsi(:)=0.0d0

      lplawv = 0
      if (ipr .ge. 50) lplawv = 1
Cdelw1       call pratfs(spid,lplawv,z,a,nr,rmax,nrmt,lmxa,pl,nsp,v,rofi,w(og),
Cdelw1      .  w(opsi))
       call pratfs ( spid , lplawv , z , a , nr , rmax , nrmt , lmxa 
     . , pl , nsp , v , rofi , rv_w_og , rv_w_opsi ) 

Cdelw1 rlse name= og old_list= og opsi 
Cdelw1 rlse name= og new_list= (None)
Cdelw1       call rlse(og)
       if (allocated(rv_w_opsi)) deallocate(rv_w_opsi)
       if (allocated(rv_w_og)) deallocate(rv_w_og)

      endif

C --- Optimise smooth-Hankel basis ---
      call dvset(rtab,1,n0,-1d0)
      call dpzero(etab,n0)
      i = 1
      if (.not. cmdopt('--noopt',7,0,strn)) then
      if (cmdopt('--norscnst',10,0,strn)) i = 0
      if (z .gt. 0) then
        call optfab(i,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,nsp,v,
     .  rofi,spid,itab,rtab,etab)
C   ... Fit value and slope of local orbitals
        call ftfalo(i,z,a,nr,rmax,nrmt,rmt,lmxa,pnu,pz,rs3,eh3,vmtz,
     .    nsp,v,rofi,spid)
      endif
      endif
      endif

C --- Print charges within/outside MT sphere ---
      if (z .gt. 0) then
      qvt = 0d0
      qct = 0d0
      qvin = 0d0
      qcin = 0d0
      do  30  isp = 1, nsp
      do  30  i = 1, nr
        r = rofi(i)
        wt0 = rofi(i+nr)
C       wt0 = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
        wt1 = wt0
        if (i .gt. nrmt) wt1 = 0d0
        if (i.eq.1 .or. i.eq.nrmt) wt1 = a*(r+b)/3d0
        if (i.eq.1 .or. i.eq.nr)   wt0 = a*(r+b)/3d0
        qvin = qvin + wt1*rho(i+(isp-1)*nr)
        qcin = qcin + wt1*rhoc(i+(isp-1)*nr)
        qvt = qvt + wt0*rho(i+(isp-1)*nr)
        qct = qct + wt0*rhoc(i+(isp-1)*nr)
   30 continue
      qtt = qvt + qct
      qtin = qvin + qcin
      if (ipr .ge. 40) write (stdo,550) qvin,qvt-qvin,qvt,qcin,qct-qcin,
     .  qct,qtin,qtt-qtin,qtt
  550 format(/' Charges:     inside',7x,'outside',7x,'sum'
     .   /' valence',3f13.6/' core   ',3f13.6/' total  ',3f13.6)

      write (stdl,710) z,rmax,qc,qct-qcin,dq,etot
  710 format('fa Z',f6.1,'   rm',f7.2,'  qc',f6.2,'  qspl',f8.5,
     .   '  dq',f8.5,'  Etot',f15.6)
      else
        qvt = 0
      endif

C --- Attach smooth Hankel tails to valence density ---
C     lplfa = nglob('lplfa')
      lplfa = 0
      if (qvt .gt. 1d-6) then
        if (lplfa .eq. 1) then
          write(stdo,344)
  344     format(/' write plot file with valence density..')
          if (is .lt. 10) write (str,'(''pl'',i1)') is
          if (is .ge. 10) write (str,'(''pl'',i2)') is
          ifi = fopna(str,-1,0)
          write (ifi,490) spid,rmt,rsmfa,nxi
  490     format('# fit to fa density: ',a/
     .       '# rmt=',f7.3,'   rsm=',f7.3,'   nxi=',i2)
          call fclose(ifi)
        endif
        call tailsm(0,nr,nrmt,nsp,a,b,rmt,rsmfa,nxi0,nxi,exi,rofi,
     .     rho,rhot,hfc,hfct)
C       call prrmsh('rho-fa',rofi,rho,nr,nr,1)
      else
        call dpzero(hfc, nxi0*nsp)
        call dpzero(hfct, nxi0*nsp)
      endif

C --- Fit analytical expression to tail of core density ---
      call fctail(nr,nrmt,a,b,rfoca,rofi,rhoc,ccof,ceh)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine freats 


      subroutine pratfs(spid,lplawv,z,a,nr,rmax,nrmt,lmaxa,pl,nsp,v,
     .  rofi,g,psi)

C- Prints out core and valence energy levels of free-atom
C ----------------------------------------------------------------------
Ci Inputs
Ci   spid  :species label
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rmax  :free-atom radius, in a.u.
Ci   nrmt  :mesh for MT radius
Ci   lmaxa :muffin-tin cutoff
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   g     :normalized wave function times r (work array)
Co Outputs
Co   psi   :normalized wave functions for each l
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n0,nr,lmaxa,nsp,nrmt,lplawv
      parameter (n0=10)
      double precision z,a,rmax,pl(n0,nsp),v(nr,nsp),rofi(nr),
     .  g(2*nr),psi(nr,0:lmaxa,nsp)
C ... Local parameters
      integer fopna,nglob,isp,l,konfig,nn,nre,i,konf,
     .  konfg(0:8),ifi,stdo,lmaxc
      double precision ev(0:20),pi,b,tol,eb1,eb2,dl,val,slo,sum,pzero,
     .  pmax,ctp,ecor,rmt
      character*1 lsym(0:n0-1), cc, str*15, spid*8
      data lsym /'s','p','d','f','g','5','6','7','8','9'/

      stdo = nglob('stdo')
      pi   = 4d0*datan(1d0)
      if (lmaxa.gt.n0-1) call rx('pratfs:  lmax too large')
      b = rmax/(dexp(a*nr-a)-1d0)
      rmt = b*(dexp(a*nrmt-a)-1d0)
      tol = 1d-8
      write(stdo,580)
  580 format(/' Free-atom wavefunctions:')

      do  80  isp = 1, nsp

C --- Valence states ---
      if (isp .eq. 1) write(stdo,401)
      if (isp .eq. 2) write(stdo,'(/'' spin 2:'')')
      eb1 = -50d0
      eb2 =  50d0
      do  20  l = 0, lmaxa
        konfig = pl(l+1,isp)
        dl = dtan(pi*(0.5d0-pl(l+1,isp)))
        nn = konfig-l-1
        ev(l) = -0.5d0
        val = rmax
        slo = dl+1
        if (rmax .gt. 9.99d0) then
          val = 1d-30
          slo = -val
        endif
        call rseq(eb1,eb2,ev(l),tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b,
     .    rofi,nr,nre)
        call gintsl(g,g,a,b,nr,rofi,sum)
        call gintsl(g,g,a,b,nrmt,rofi,pmax)
        sum = sum - pmax
        call ppratf(ev(l),z,nr,nre,rofi,a,b,v(1,isp),g,pzero,pmax,ctp)
        cc = ' '
        if (dabs(ctp-rmax) .lt. 1d-3) cc = '*'
        write(stdo,400) konfig,lsym(l),ev(l),pzero,pmax,ctp,cc,sum
  400   format(i4,a1,f14.5,2x,3f12.3,a,f12.6)
  401   format(' valence:',6x,'eval',7x,'node at',6x,'max at',7x,
     .    'c.t.p.   rho(r>rmt)')

C   ... Copy valence wavefunction to psi
        do  24  i = 1, nr
   24   psi(i,l,isp) = g(i)
   20 continue

C --- Core states ---
      write(stdo,403)
      eb1 = -2.5d0*z*z-5d0
      eb2 = 50d0
      call config(pl,lmaxa,z,konfg,lmaxc)

      do  40  konf = 1, 8
      do  40  l = 0, min(konf-1,lmaxc)
      konfig = konfg(l)
      if (konf .ge. konfig) goto 40
      nn = konf-l-1
      ecor = -50d0
      val = 1d-30
      slo = -val
      call rseq(eb1,eb2,ecor,tol,z,l,nn,val,slo,v(1,isp),g,sum,a,b,rofi,
     .  nr,nre)
      call gintsl(g,g,a,b,nr,rofi,sum)
      call gintsl(g,g,a,b,nrmt,rofi,pmax)
      sum = sum - pmax
      call ppratf(ecor,z,nr,nre,rofi,a,b,v(1,isp),g,pzero,pmax,ctp)
      write(stdo,400) konf,lsym(l),ecor,pzero,pmax,ctp,' ',sum
  403 format(/' core:        ecore',7x,'node at',6x,'max at',
     .   7x,'c.t.p.   rho(r>rmt)')
  40  continue

  80  continue

C --- Write file with valence wavefunctions
      if (lplawv .eq. 1) then
        write (str,'(''wf_'',a)') spid
        write (stdo,344) str
  344   format(/' Write valence wavefunctions to plot file: ',a)
        ifi = fopna(str,-1,0)
        write (ifi,490) spid,rmax,rmt,nr,lmaxa,nr,1+nsp*(lmaxa+1)
  490   format('# Free-atom wavefunctions (divided by r) for species ',
     .     a/'# rmax=',f7.3,'   rmt=',f7.3,'   nr=',i5,'   lmax=',i3/
     .    '% rows ',i5,' cols ',i3)
        do  50  i=1,nr
          write (ifi,495) rofi(i),((psi(i,l,isp),l=0,lmaxa),isp=1,nsp)
  495     format(f9.5,1p,16d14.5)
   50   continue
        call fclr(str,ifi)
      endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pratfs 


      subroutine ppratf(e,z,nr,nre,rofi,a,b,v,g,pzero,pmax,ctp)

C- Find outermost node and maximum of wavefct
C ----------------------------------------------------------------------
Ci Inputs
Ci   e     :wave function eigenvalue
Ci   z     :nuclear charge
Ci   nr    :number of radial mesh points
Ci   nre   :last point for which wf is calculated
Ci   rofi  :radial mesh points
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   v     :spherical potential (atomsr.f)
Ci   g     :normalized wave function times r
Co Outputs
Co   pzero :outermost node
Co   pmax  :outermost maximum
Co   ctp   :classical turning point
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nre
      double precision a,b,ctp,e,pmax,pzero,z,rofi(nr),v(nr),g(nr)
C ... Local parameters
      integer i,ir
      double precision g1,g2,rho1,rho2,rho3,x

C ... Find the classical turning point
      do  20 i = nr-1, 5, -1
        ir = i
        if (e .gt. v(i)-2d0*z/rofi(i)) goto 21
   20 continue
   21 g1 = e-v(ir) + 2d0*z/rofi(ir)
      g2 = e-v(ir+1) + 2d0*z/rofi(ir+1)
      ctp = rofi(nr)
      if (g1*g2 .lt. 0d0) ctp = (rofi(ir)*g2-rofi(ir+1)*g1)/(g2-g1)

C ... Find the outermost node
      do  10  i = nre-1, 5, -1
        ir = i
        if (g(i)*g(i+1) .lt. 0d0) goto 11
   10 continue
   11 continue
      pzero = 0d0
      g1 = g(ir)
      g2 = g(ir+1)
      if (ir .gt. 5) pzero = (rofi(ir)*g2-rofi(ir+1)*g1)/(g2-g1)

C ... Find the outermost maximum
      do  30  i = nre-2, 5, -1
        ir = i
        rho1 = g(i)*g(i)
        rho2 = g(i+1)*g(i+1)
        rho3 = g(i+2)*g(i+2)
        if (rho1 .lt. rho2) goto 31
   30 continue
   31 pmax = 0
      if (ir .gt. 5) then
        x = -0.5d0*(rho3-rho1)/(rho1+rho3-2*rho2)
        pmax = b*(dexp(a*(ir+x))-1d0)
      endif
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine ppratf 


      subroutine optfab(isw,z,a,nr,rmax,nrmt,rmt,lmxa,pl,ql,nsp,v,rofi,
     .  spid,itab,rtab,etab)

C- Optimise a minimal smooth-Hankel basis for the free atom.
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :1 constrain rsm to be <= rmt
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rmax  :muffin-tin radius, in a.u.
Ci   nrmt  :number of points between 0..rmt
Ci   rmt   :muffin-tin radius, in a.u.
Ci   lmxa  :muffin-tin l-cutoff
Ci   pl    :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   ql    :sphere moments
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   spid  :species label
Co Outputs
Co   itab  :itab(l+1)=1  optimized wave function was found for this l
Co   rtab  :smoothing radius for optimized wave function
Co   etab  :energy for optimized wave function
Cl Local variables
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmxa,nr,nrmt,nsp,n0,isw
      double precision a,rmax,rmt,z
      parameter (n0=10)
      double precision rofi(1),v(nr,nsp),pl(n0,nsp),ql(3,n0,nsp)
      character spid*8
C ... Local parameters
      logical cmdopt
      character strn*80
      integer itab(n0,2)
      double precision rtab(n0,2),etab(n0,2)
Cdelw1       integer ipr,iprint,irep,isp,istife,istifr,jpr,konfig,l,
Cdelw1      .  lgunit,lplawv,lrel,nglob,nn,nrep,og,ogp,oh,opsi,stdo,stdl
       integer:: ipr , iprint , irep , isp , istife , istifr , jpr , 
     . konfig , l , lgunit , lplawv , lrel , nglob , nn , nrep , stdo 
     . , stdl 
       real(8) ,allocatable :: rv_w_og(:)
       real(8) ,allocatable :: rv_w_ogp(:)
       real(8) ,allocatable :: rv_w_oh(:)
       real(8) ,allocatable :: rv_w_opsi(:)

      double precision b,deh,deh0,dphi,dphip,drsm,drsm0,e1,e2,e3,eadd,
     .  eaddx,eh,elim1,elim2,enew,enu,eval,p,phi,phip,pnu,qvl,radd,
     .  raddx,rlim1,rlim2,rnew,rsm,stife,stifr,sume1,sume2,qrmt
C ... Heap
      integer w(1)
      common /w/ w

      ipr = iprint()
      stdo = lgunit(1)
      stdl = lgunit(2)
      lrel = nglob('lrel')
      if (z .lt. 0.99d0) lrel = 0
      if (lmxa .gt. 8) call rx('optfab:  lmax too large')
      b = rmax/(dexp(a*nr-a)-1d0)

      do  80  isp = 1, nsp
      if (ipr .ge. 20) call awrit2('%?#n>=30#%N##'//
     .  ' Optimise free-atom basis for species '//spid//'%a, rmt=%;7g',
     .  ' ',80,stdo,ipr,rmt)
Cdelw1       call defrr (oh,    nr)
       allocate(rv_w_oh(nr))
       if (nr<0) rv_w_oh(:)=0.0d0

Cdelw1       call defrr (og,    2*nr)
       allocate(rv_w_og(2*nr))
       if (2*nr<0) rv_w_og(:)=0.0d0

Cdelw1       call defrr (ogp,   2*nr*4)
       allocate(rv_w_ogp(2*nr*4))
       if (2*nr*4<0) rv_w_ogp(:)=0.0d0


C --- Parameters for minimisation ---
      drsm0 = 0.1d0
c takao makes "safer setting"
ctakao
c      rlim1 = 0.3d0
      rlim1 = 0.9d0
ctakao
c      rlim2 = 2*rmt
      rlim2 = rmt+1d-6
 
      raddx = 0.2d0

      deh0  = 0.05d0
ctakao
c      elim1 = -5.0d0
      elim1 = -2.5d0
      elim2 = -0.10d0
C     elim2 = -0.20d0
      eaddx = 0.099d0
      jpr=0
      if (ipr .ge. 50) jpr=1

C --- Loop over bound valence states ---
      if (ipr .ge. 20) write (stdo,261)
      sume1 = 0d0
      sume2 = 0d0
      do  10  l = 0, lmxa
        itab(l+1,isp) = 0
        konfig = pl(l+1,isp)
        nn = konfig-l-1
        qvl = ql(1,l+1,isp)
C   ... Get exact fa wavefunction, eigval, pnu at rmt
Cdelw1         call popta3(0,l,z,nn,rmt,nr,nrmt,rofi,v(1,isp),a,b,
Cdelw1      .    eval,pnu,w(og))
         call popta3 ( 0 , l , z , nn , rmt , nr , nrmt , rofi , v ( 1 
     .   , isp ) , a , b , eval , pnu , rv_w_og ) 

        if (eval .gt. 0d0) goto 10
        sume1 = sume1 + qvl*eval
C   ... Potential parameters at MT sphere
Cdelw1         call popta4(l,z,rmt,nrmt,rofi,v(1,isp),w(og),w(ogp),
Cdelw1      .    a,b,pnu,enu,p,phi,dphi,phip,dphip)
         call popta4 ( l , z , rmt , nrmt , rofi , v ( 1 , isp ) , rv_w_og 
     .   , rv_w_ogp , a , b , pnu , enu , p , phi , dphi , phip , dphip 
     .   ) 

        rsm = rmt
        eh = -1
        if (jpr .gt. 0) write (stdo,340)
  340   format('  L   parin    aux      E1       E2       E3',
     .     '       stiff    Eout     parout')
        do  12  irep = 1, 50
          nrep = irep
C     ... Get center energy
Cdelw1           call popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e2,qrmt)
           call popta1 ( rsm , eh , l , z , rmt , nr , nrmt , rofi , rv_w_oh 
     .     , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip , dphip 
     .     , e2 , qrmt ) 

C     ... Vary rsm
          drsm = drsm0
Cdelw1           call popta1(rsm+drsm,eh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e3,qrmt)
           call popta1 ( rsm + drsm , eh , l , z , rmt , nr , nrmt , rofi 
     .     , rv_w_oh , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip 
     .     , dphip , e3 , qrmt ) 

Cdelw1           call popta1(rsm-drsm,eh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e1,qrmt)
           call popta1 ( rsm - drsm , eh , l , z , rmt , nr , nrmt , rofi 
     .     , rv_w_oh , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip 
     .     , dphip , e1 , qrmt ) 

          call popta2(l,rsm,eh,drsm,e1,e2,e3,rlim1,rlim2,raddx,rnew,
     .       stifr,jpr)
C     ... Vary eh
          deh = deh0
c         if (eh+deh.gt.-0.01d0) deh=-eh-0.01d0
Cdelw1           call popta1(rsm,eh+deh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e3,qrmt)
           call popta1 ( rsm , eh + deh , l , z , rmt , nr , nrmt , rofi 
     .     , rv_w_oh , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip 
     .     , dphip , e3 , qrmt ) 

Cdelw1           call popta1(rsm,eh-deh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e1,qrmt)
           call popta1 ( rsm , eh - deh , l , z , rmt , nr , nrmt , rofi 
     .     , rv_w_oh , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip 
     .     , dphip , e1 , qrmt ) 

          call popta2(l,eh,rsm,deh,e1,e2,e3,elim1,elim2,eaddx,enew,
     .       stife,jpr)

          radd = rnew-rsm
          eadd = enew-eh
          rsm = rnew
          eh = enew
          if (dabs(radd).lt.5d-3 .and. dabs(eadd).lt. 5d-3) goto 90
   12   continue
   90   continue

C   ... End of iteration loop

        sume2 = sume2 + qvl*e2
        if (ipr .ge. 20)
     .  write (stdo,260) l,nrep,rsm,eh,stifr,stife,e2,eval,pnu,qvl
  260   format(i2,i4,2f8.3,1x,2f9.1,1x,2f10.5,f8.2,f7.2)
  261   format(' l  it    Rsm      Eh     stiffR   stiffE',
     .     '      Eval      Exact     Pnu    Ql')
        istifr = stifr+0.5d0
        istife = stife+0.5d0
        write (stdl,710) l,nrep,rsm,eh,istifr,istife,e2,eval,pnu,qvl
  710   format('fa op',i2,i4,2f7.3,'  stf',2i6,'  ev',2f9.5,
     .     '  pq',2f6.2)

C   ... Possibly constrain rsm
ctakao
        if (mod(isw,10) .eq. 1 .and. rsm .gt. rmt*2d0/3d0) then
        if (ipr .ge. 20)
     .  write(stdo,
c takao '('' ... rsm exceeded rmt .. repeat with rsm=rmt'')')
     .  '('' ... rsm exceeded rmt*2/3 .. repeat with rsm=rmt'')')
c takao
c        rsm = rmt
        rsm = rmt*2d0/3d0
        sume2 = sume2 - qvl*e2

        do  112  irep = 1,50
          nrep = irep
C     ... Get center energy
Cdelw1           call popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e2,qrmt)
           call popta1 ( rsm , eh , l , z , rmt , nr , nrmt , rofi , rv_w_oh 
     .     , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip , dphip 
     .     , e2 , qrmt ) 

C     ... Vary eh
          deh = deh0
c         if (eh+deh.gt.-0.01d0) deh=-eh-0.01d0
Cdelw1           call popta1(rsm,eh+deh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e3,qrmt)
           call popta1 ( rsm , eh + deh , l , z , rmt , nr , nrmt , rofi 
     .     , rv_w_oh , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip 
     .     , dphip , e3 , qrmt ) 

Cdelw1           call popta1(rsm,eh-deh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .       v(1,isp),a,b,enu,p,phi,dphi,phip,dphip,e1,qrmt)
           call popta1 ( rsm , eh - deh , l , z , rmt , nr , nrmt , rofi 
     .     , rv_w_oh , v ( 1 , isp ) , a , b , enu , p , phi , dphi , phip 
     .     , dphip , e1 , qrmt ) 

          call popta2(l,eh,rsm,deh,e1,e2,e3,elim1,elim2,eaddx,enew,
     .       stife,jpr)

          eadd = enew-eh
          eh = enew
          if (dabs(eadd).lt. 5d-3) goto 190
  112   continue
  190   continue
C   ... End of iteration loop

        sume2 = sume2 + qvl*e2
        if (ipr .ge. 20)
     .  write (stdo,260) l,nrep,rsm,eh,stifr,stife,e2,eval,pnu,qvl
        istife = stife+0.5d0
        write (stdl,710) l,nrep,rsm,eh,istifr,istife,e2,eval,pnu,qvl
        endif

        itab(l+1,isp) = 1
        rtab(l+1,isp) = rsm
        etab(l+1,isp) = eh

   10 continue

      if (ipr .ge. 20) write (stdo,320) sume1,sume2,sume2-sume1
  320 format(' eigenvalue sum:  exact',f10.5,'    opt basis',f10.5,
     .   '    error',f8.5)
      write (stdl,720) sume1,sume2,sume2-sume1
  720 format('fa op sumev',f11.5,'   opt basis',f11.5,'   err',f9.5)

   80 continue

C --- Make plot file ---
C     lplawv=nglob('lplawv')
      lplawv = 0
      if (cmdopt('--plotwf',8,0,strn)) lplawv = 1
      if (lplawv .eq. 1) then
        if (nsp .eq. 2) call rx('optfab is not spinpol yet')
Cdelw1         call defrr(opsi,  nr*lmxa)
         allocate(rv_w_opsi(nr*lmxa))
         if (nr*lmxa<0) rv_w_opsi(:)=0.0d0

Cdelw1         call popta5(lmxa,rtab,etab,itab,z,pl,rmax,rmt,nr,nrmt,
Cdelw1      .     rofi,w(opsi),v,w(og),a,b,spid)
         call popta5 ( lmxa , rtab , etab , itab , z , pl , rmax , rmt 
     .   , nr , nrmt , rofi , rv_w_opsi , v , rv_w_og , a , b , spid ) 

Cdelw1 rlse name= opsi old_list= oh og ogp opsi 
Cdelw1 rlse name= opsi new_list= oh og ogp 
Cdelw1         call rlse(opsi)
         if (allocated(rv_w_opsi)) deallocate(rv_w_opsi)

      endif

Cdelw1 rlse name= oh old_list= oh og ogp 
Cdelw1 rlse name= oh new_list= (None)
Cdelw1       call rlse (oh)
       if (allocated(rv_w_ogp)) deallocate(rv_w_ogp)
       if (allocated(rv_w_og)) deallocate(rv_w_og)
       if (allocated(rv_w_oh)) deallocate(rv_w_oh)


Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine optfab 


      subroutine ftfalo(icst,z,a,nr,rmax,nrmt,rmt,lmxa,pnu,pz,rs3,eh3,
     .  vmtz,nsp,v,rofi,spid)

C- Fit value and slope of local orbitals to smoothed Hankel
C ----------------------------------------------------------------------
Ci Inputs
Ci   icst  :1 constrain rsm to be <= rmt
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(ir) = b [e^(a(ir-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rmax  :muffin-tin radius, in a.u.
Ci   nrmt  :number of points between 0..rmt
Ci   rmt   :muffin-tin radius, in a.u.
Ci   lmxa  :muffin-tin l-cutoff
Ci   pl    :boundary conditions for valence wavefunctions.
Ci   pz    :boundary conditions for local orbital. pz=0 -> no loc. orb.
Ci         :10s digit controls how local orbital included in hamiltonian
Ci         :10s digit nonzero -> smooth Hankel tail is attached.
Ci   rs3   :minimum allowed smoothing radius in attaching Hankel tails
Ci         :to local orbitals
Ci   eh3   :Hankel energy when attaching Hankel tails to high-lying
Ci         :local orbitals
Ci   vmtz  :parameter used in attaching Hankel tails to local orbitals
Ci         :It is used as a constant shift to Hankel energies for the 
Ci         :fitting of local orbitals to Hankel tails. Thus vmtz
Ci         :is an estimate for the potential at the MT radius.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   spid  :species label
Co Outputs
Cl Local variables
Cr Remarks
Cu Updates
Cu   16 Jun 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmxa,nr,nrmt,nsp,n0,icst
      double precision a,rmax,rmt,z,rs3,eh3,vmtz
      parameter (n0=10)
      double precision rofi(1),v(nr,nsp),pz(n0,nsp),pnu(n0,nsp)
C ... Local parameters
      logical cmdopt
      character spid*8, strn*80, flg(2)*1
Cdelw1       integer ipr,iprint,i,konfig,l,info,lgunit,nn,og,ogp,oh,stdo,
Cdelw1      .  lplawv,loclo,nfit,isw
       integer:: ipr , iprint , i , konfig , l , info , lgunit , nn 
     . , stdo , lplawv , loclo , nfit , isw 
       real(8) ,allocatable :: rv_w_og(:)
       real(8) ,allocatable :: rv_w_ogp(:)
       real(8) ,allocatable :: rv_w_oh(:)

Cdelw1       integer opsi
       real(8) ,allocatable :: rv_w_opsi(:)

      double precision b,dasum,dphi,dphip,e2,eh,eval,p,phi,phip,
     .  pnul,rsm,rsmin,rsmax,ekin
C     emin and emax are the maximum allowed ranges in Hankel energies
C     for the fitting of local orbitals to Hankel tails.
      double precision emin,emax,tphi
C     For plotting wave functions
      integer itab(n0,2)
      double precision rtab(n0,2),etab(n0,2),pl(n0,nsp),qrmt
C ... Heap
      integer w(1)
      common /w/ w
      data flg/'*',' '/

C     return

      ipr = iprint()
      stdo = lgunit(1)
C     stdl = lgunit(2)
      if (lmxa .gt. 8) call rx('ftfalo:  lmax too large')
      b = rmax/(dexp(a*nr-a)-1d0)
      nfit = 0

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (dasum(lmxa+1,pz,1) .eq. 0) return

      do  80  i = 1, nsp
Cdelw1       call defrr (oh,    nr)
       allocate(rv_w_oh(nr))
       if (nr<0) rv_w_oh(:)=0.0d0

Cdelw1       call defrr (og,    2*nr)
       allocate(rv_w_og(2*nr))
       if (2*nr<0) rv_w_og(:)=0.0d0

Cdelw1       call defrr (ogp,   2*nr*4)
       allocate(rv_w_ogp(2*nr*4))
       if (2*nr*4<0) rv_w_ogp(:)=0.0d0


C --- Loop over local orbitals ---
C      sume1 = 0d0
C      sume2 = 0d0
      do  10  l = 0, lmxa

        itab(l+1,i) = 0
        pnul = pnu(l+1,i)
        pl(l+1,i) = pnu(l+1,i)
        konfig = mod(pz(l+1,1),10d0)

C       Skip all but local orbitals with tails attached
        if (pz(l+1,1) .lt. 10) goto 10

C       Case local orbital deeper than valence
        if (int(pnul-1) .eq. int(mod(pz(l+1,1),10d0))) then
          loclo = 1
C         Not needed, actually, since overwritten by popta3
C         pnul = mod(pz(l+1,1),10d0)

C       Case local orbital higher than the valence state
        elseif (int(pnul+1) .eq. int(mod(pz(l+1,1),10d0))) then
          pnul = mod(pz(l+1,1),10d0)
          loclo = 0

C       Local orbital neither one: error
        else
          call fexit3(-1,111,' Exit -1 freeat, l=%i:  sc '//
     .      'PZ=%d incompatible with valence P=%;3d',l,pz(l+1,1),pnul)
        endif

C       Skip high-lying local orbitals unless specifically sought
        if (loclo .eq. 0 .and. .not. cmdopt('--getallloc',11,0,strn))
     .    goto 10

        nfit = nfit + 1
        if (nfit .eq. 1) then
          call info2(20,1,0,
     .    ' Fit local orbitals to sm hankels, species '//spid//
     .    '%a, rmt=%;7g',rmt,0)
          if (ipr .ge. 20) write (stdo,261)
        endif

C   ... Get exact fa wavefunction, eigval, pnu_l at rmt
        if (loclo .eq. 1) then
          nn = konfig-l-1
Cdelw1           call popta3(0,l,z,nn,rmt,nr,nrmt,rofi,v(1,i),a,b,eval,
Cdelw1      .      pnul,w(og))
           call popta3 ( 0 , l , z , nn , rmt , nr , nrmt , rofi , v ( 1 
     .     , i ) , a , b , eval , pnul , rv_w_og ) 

C       Finish if in future, need w.f. at r>rmt
C        else
C          call popta3(1,l,z,nn,rmt,nr,nrmt,rofi,v(1,i),a,b,eval,
C     .      pnul,w(og))
        endif
        pl(l+1,i) = pnul

C   ... Potential parameters at MT sphere
Cdelw1         call popta4(l,z,rmt,nrmt,rofi,v(1,i),w(og),w(ogp),
Cdelw1      .    a,b,pnul,eval,p,phi,dphi,phip,dphip)
         call popta4 ( l , z , rmt , nrmt , rofi , v ( 1 , i ) , rv_w_og 
     .   , rv_w_ogp , a , b , pnul , eval , p , phi , dphi , phip , dphip 
     .   ) 


C   ... Set conditions on envelope functions ... For now
        rsmin = rs3
        rsmax = 5
        if (icst .eq. 1) rsmax = rmt
C       Use r->infty value for energy
        eh = min(-.02d0,eval)

C   ... Match Hankel to phi,dphi
C        rsm = rsmin
C        emax = -.02d0
C        emin = -5d0
C        call mtchre(100,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,phi,
C     .    dphi,rsm,eh,ekin,info)

C   ... Match slope and K.E. of Hankel to phi,dphi
        tphi = eval - (v(nrmt,i)-2*z/rmt)
        rsm = 0
        eh = min(eval-vmtz,-.02d0)
        emax = -.02d0
        emin = -10d0
C       if (ipr .ge. 20) call pshpr(max(ipr,50))
        call mtchre(003,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,tphi,
     .    dphi,rsm,eh,ekin,info)
C       if (ipr .ge. 20) call poppr
C       Match failed ... turn up verbosity and repeat for info
        if (info .eq. -1) then
          call info2(0,2,1,
     .      ' *** ftfalo (fatal) cannot fit smooth Hankel to w.f.'//
     .      ' class '//spid//
     .      '%N ... possibly reduce RS3 (current value = %,1d)',rs3,0)
          call pshpr(max(ipr,110))
          call mtchr2(1,l,emin,emax,(emin+emax)/2,
     .      rmt,phi,dphi,rsmin,eh,ekin,i)
          call poppr
C         call pshpr(max(ipr,110))
C         call mtchre(103,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,tphi,
C    .      dphi,rsm,eh,ekin,info)
          call fexit2(-1,111,
     .      ' Exit -1 : ftfalo : failed to match log der=%,1;3d'//
     .      ' to envelope, l=%i',dphi/phi,l)
        endif

C  ... Get energy of this wave function
Cdelw1         call popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,w(oh),
Cdelw1      .    v(1,i),a,b,eval,p,phi,dphi,phip,dphip,e2,qrmt)
         call popta1 ( rsm , eh , l , z , rmt , nr , nrmt , rofi , rv_w_oh 
     .   , v ( 1 , i ) , a , b , eval , p , phi , dphi , phip , dphip 
     .   , e2 , qrmt ) 


        if (ipr .ge. 20)
     .    write (stdo,260) l,rsm,eh,qrmt,e2,eval,pnul,tphi,ekin,
     .    flg(2-isw(dabs(ekin-tphi).gt.1d-5))

  260   format(i2,2f8.3,3f10.5,f9.3,2f10.5,a1,f10.5)
  261   format(' l    Rsm     Eh     Q(r>rmt)   Eval',
     .    '      Exact      Pnu     K.E.    fit K.E.')

        itab(l+1,i) = 1
        rtab(l+1,i) = rsm
        etab(l+1,i) = eh

   10 continue
   80 continue

C --- Make plot file ---
C     lplawv=nglob('lplawv')
      lplawv = 0
      if (cmdopt('--plotwf',8,0,strn)) lplawv = 1
      if (lplawv .eq. 1) then
        if (nsp .eq. 2) call rx('optfab is not spinpol yet')
Cdelw1         call defrr(opsi,  nr*lmxa)
         allocate(rv_w_opsi(nr*lmxa))
         if (nr*lmxa<0) rv_w_opsi(:)=0.0d0

Cdelw1         call popta5(lmxa,rtab,etab,itab,z,pl,rmax,rmt,nr,nrmt,
Cdelw1      .     rofi,w(opsi),v,w(og),a,b,spid)
         call popta5 ( lmxa , rtab , etab , itab , z , pl , rmax , rmt 
     .   , nr , nrmt , rofi , rv_w_opsi , v , rv_w_og , a , b , spid ) 

Cdelw1 rlse name= opsi old_list= oh og ogp opsi 
Cdelw1 rlse name= opsi new_list= oh og ogp 
Cdelw1         call rlse(opsi)
         if (allocated(rv_w_opsi)) deallocate(rv_w_opsi)

      endif

Cdelw1 rlse name= oh old_list= oh og ogp 
Cdelw1 rlse name= oh new_list= (None)
Cdelw1       call rlse (oh)
       if (allocated(rv_w_ogp)) deallocate(rv_w_ogp)
       if (allocated(rv_w_og)) deallocate(rv_w_og)
       if (allocated(rv_w_oh)) deallocate(rv_w_oh)


Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine ftfalo 


      subroutine popta1(rsm,eh,l,z,rmt,nr,nrmt,rofi,h,v,a,b,enu,p,
     .  phi,dphi,phip,dphip,eval,qrmt)

C- Calculate expectation value for smooth Hankel
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsm   :smoothing radius of basis function
Ci   eh    :energy of basis function
Ci   l     :l quantum number
Ci   z     :nuclear charge
Ci   rmt   :muffin-tin radius
Ci   nr    :number of radial mesh points
Ci   nrmt  :number points to muffin-tin radius
Ci   rofi  :radial mesh points
Ci   h     :work array
Ci   v     :spherical potential (atomsr.f)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   enu   :enu's for making charge density
Ci   p     :<gp**2> (potential parameter)
Co   phi   :wave function at rmt
Co   dphi  :radial derivative of of phi at rmt
Co   phip  :energy derivative of phi
Co   dphip :radial derivative of dphi
Co Outputs
Co   eval  :expectation value
Co   qrmt  :fraction of (wave function)^2 for r>rmt
Cr Remarks
Cu Updates
Cu   24 Sep 04 return qrmt
Cu   16 Jun 04 Adapted to new hansmd, mtchae
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,nr,nrmt
      double precision a,b,dphi,dphip,eh,enu,eval,p,phi,phip,rmt,rsm,z,
     .  rofi(nr),h(nr),v(nr),qrmt
C ... Local parameters
      integer i
      double precision alfa,beta,det,drdi,hlap,
     .  hum,hum1,hum2,r,sum,sum1,sum2,tum2,vum2,wt
C     double precision xi(0:20)

      double precision hs(0:l),dhs(0:l),ddhs(0:l)

C     pi = 4d0*datan(1d0)
C     asm = 1d0/rsm
C     lp1 = l+1

C ... Integrals over smooth Hankel on mesh
C     gfac = (asm*asm/pi)**1.5d0 * dexp(eh*rsm*rsm/4d0)
C     ta2 = 2d0*asm*asm
      tum2 = 0d0
      sum2 = 0d0
      vum2 = 0d0

      do  10  i = nrmt, nr
        r = rofi(i)

C   ... Make r*h and r Laplacian h, including L^2
        call hansmd(2,r,eh,rsm,l,hs,dhs,ddhs,det,det,det)
        h(i) = hs(l)*r
        hlap = ddhs(l)*r
CC      Old : r*h and r Laplacian h, including L^2
CC      h = r*radial part of sm. Hankel
C       call hansmr(r,eh,asm,xi,l)
C       h(i) = xi(l)*(r**lp1)
CC      radial part of Gaussian
C       gl = gfac * dexp(-asm*asm*r*r) * ta2**l * (r**lp1)
CC      r * (nabla_r - l(l+1)/r^2) h_l
C       hlap = -4d0*pi*gl - eh*h(i)

C  ...  Accumulate <h h>, <h v h>, <h -nabla h>
        wt = 2*(mod(i+1,2)+1)/3d0
        if (i.eq.nrmt .or. i.eq.nr) wt = 1d0/3d0
        drdi = a*(r+b)
        sum2 = sum2 + wt*drdi*h(i)*h(i)
        vum2 = vum2 + wt*drdi*h(i)*h(i)*(v(i)-2d0*z/r)
        tum2 = tum2 + wt*drdi*h(i)*(-hlap)
   10 continue
      hum2 = tum2+vum2

C --- BC's: match phi,phidot to envelope at RMT ---
      call mtchae(0,rsm,eh,l,rmt,phi,dphi,phip,dphip,alfa,beta)
CC    OLD matching
CC    Match value, slope fl,dfl to linear combination of phi,phidot
C     call hansmr(rmt,eh,asm,xi,l+1)
CC    Value and radial derivative of h (JMP 39, 3393, Eq. 4.7)
C     fl = xi(l)*rmt**l
C     flp1 = xi(l+1)*rmt**(l+1)
C     dfl = l*fl/rmt-flp1
CC    Match fl,dfl to linear combination of phi,phidot
CC    Use  phi=phi(R); phip=phidot(R) dphi=phi'(R); dphip=phidot'(R)
CC    (phi  phip ) (alpha)   (fl )    (alpha)    1  (dphip -phip) (fl )
CC    (          ) (     ) = (   ) -> (     ) = --- (           ) (   )
CC    (dphi dphip) (beta )   (dfl)    (beta )   det (-dphi  phi ) (dfl)
C     det = phi*dphip-dphi*phip
C     alfa = (fl*dphip-dfl*phip)/det
C     beta = (dfl*phi-fl*dphi)/det

C     O = alpha^2 <phi | phi> + beta^2 <phidot | phidot>
      sum1 = alfa*alfa + beta*beta*p
      hum1 = alfa*alfa*enu + alfa*beta + beta*beta*enu*p

      sum = sum1+sum2
      hum = hum1+hum2
      eval = hum/sum

      qrmt = sum2/sum
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine popta1 


      subroutine popta2(l,x0,y0,dx,e1,e2,e3,xmin,xmax,xshx,xnew,stiff,
     .  jpr)

C- Find minimum from three values
C ----------------------------------------------------------------------
Ci Inputs
Ci   l     :angular momentum
Ci   x0    :starting value
Ci   y0    :used for printout
Ci   dx    :excursion in x for numerical differentiation
Ci   e1    :function value at x0-dx
Ci   e2    :function value at x0
Ci   e3    :function value at x0+dx
Ci   xmin  :boundary: estimated minimum must be >= xmin
Ci   xmax  :boundary: estimated minimum must be <= xmax
Ci   xshx  :maximum step size
Ci   jpr   :printout verbosity
Co Outputs
Co   xnew  :new estimate for the minimum
Co   stiff :estimated curvature
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer jpr,l
      double precision dx,e1,e2,e3,stiff,x0,xmax,xmin,xnew,xshx,y0
C ... Local parameters
      integer ie0,lgunit,stdo
      double precision a,aa,b,c,ee1,ee2,ee3,een,enew,xadd

      stdo = lgunit(1)
      c = e2
      b = (e3-e1)/(2*dx)
      a = (e1+e3-2*e2)/(2*dx*dx)
      if (a .le. 0d0) then
        xadd = -xshx
        enew = e1
        if (e3 .lt. e1) xadd = xshx
        if (e3 .lt. e1) enew = e3
      else
        xadd = -b/(2*a)
        enew = a*xadd*xadd + b*xadd + c
      endif
      aa = 2*1d3*a

      if (xadd .gt. xshx)  xadd = xshx
      if (xadd .lt. -xshx) xadd = -xshx
      xnew = x0+xadd
      if (xnew .gt. xmax) xnew = xmax
      if (xnew .lt. xmin) xnew = xmin

      ie0 = e2
      ee1 = 1d3*(e1-ie0)
      ee2 = 1d3*(e2-ie0)
      ee3 = 1d3*(e3-ie0)
      een = 1d3*(enew-ie0)
      stiff = aa

      if (jpr.gt.0) write (stdo,810)l,x0,y0,ee1,ee2,ee3,aa,een,xnew
  810 format(i3,f8.3,f8.3,f10.3,2f9.3,f9.1,f10.3,f8.3,a)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine popta2 

      subroutine popta3(mode,l,z,nn,rmt,nr,nrmt,rofi,v,a,b,evl,pnu,g)

C- Get exact fa wavefunction, eigval, pnu at Rmt.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 boundary condition is val,slo = 0 at nr
Ci         :1 boundary condition is that w.f. satisfy pnu at nrmt
Ci         :  (under development)
Ci   l     :angular momentum
Ci   z     :nuclear charge
Ci   rmax  :sphere radius
Ci   rmt   :muffin-tin radius, in a.u.
Ci   nr    :number of radial mesh points
Ci   nrmt  :number of radial mesh points to rmt
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio Inputs/Outputs
Cio  nn    :number of nodes (input mode 0; output mode 1)
Cio  pnu   :boundary condition at rmt (output mode 0; input mode 1)
Co Outputs
Co   g     :normalized wave function times r
Co   evl   :eigenvalue
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,l,nn,nr,nrmt
      double precision a,b,evl,pnu,rmt,z,rofi(nr),v(nr),g(nr*2)
C ... Local parameters
      integer lgunit,nre,stdo,konfig,nri,nn2
      double precision d0l,p0l,dphi,drdi,du,eb1,eb2,g1,g2,g3,g4,g5,pi,
     .  slo,slou,sum,tol,val,valu,dnu

      stdo = lgunit(1)
      pi = 4d0*datan(1d0)

      eb1 = -30
      eb2 = 20
      tol = 1d-10
      val = 1d-30
      slo = -val
      evl = -0.5d0
      nri = nr
      if (mode .eq. 1) then
        konfig = pnu
        nn = konfig-l-1
        dnu = dtan(pi*(0.5d0-pnu))
        val = rmt
        slo = dnu+1d0
        nri = nrmt
      endif
      call rseq(eb1,eb2,evl,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nri,nre)
      if (mode .eq. 1) then
C       integration becomes rather strange for r>>rmt.
C       Need to truncate radius.
        call rsq1(nri,evl,l,z,v,nr,g,val,slo,nn2,a,b,rofi,nr)
C       call prrmsh('g',rofi,g,nr,nr,1)
        call rx('not finished mode 1')
      endif
      g1 = g(nrmt-2)
      g2 = g(nrmt-1)
      g3 = g(nrmt)
      g4 = g(nrmt+1)
      g5 = g(nrmt+2)
      drdi = a*(rmt+b)
      valu = g3
      slou = (-2*g5+16*g4-16*g2+2*g1)/(24d0*drdi)
      du   = rmt*slou/valu
      dphi = du-1
      pnu  = nn+l+1 + (0.5d0-datan(dphi)/pi)

C ... Don't set too low..
      d0l = l
      p0l = nn+l+1 + 0.5d0-datan(d0l)/pi
      p0l = nn+l+1 + 0.1d0
      if (pnu .lt. p0l) then
        write (stdo,145) l,pnu,p0l
  145   format(' l=',i1,'  increase Pnu=',f8.3,'  to ',f8.3)
        pnu = p0l
      endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine popta3 

      subroutine popta4(l,z,rmt,nrmt,rofi,v,g,gp,a,b,pnu,enu,p,phi,dphi,
     .  phip,dphip)

C- Potential parameters at MT sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   l     :angular momentum
Ci   z     :nuclear charge
Ci   rmt   :muffin-tin radius, in a.u.
Ci   nrmt  :number of radial mesh points to rmt
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   g     :normalized wave function times r
Ci   gp    :energy derivative(s) of g
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   enu   :enu's for making charge density
Co Outputs
Co   phi   :wave function at rmt
Co   dphi  :radial derivative of of phi at rmt
Co   phip  :energy derivative of phi
Co   dphip :radial derivative of dphi
Co   p     :<gp**2> (potential parameter)
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,nrmt
      double precision a,b,dphi,dphip,enu,p,phi,phip,pnu,rmt,z
      double precision rofi(nrmt),v(nrmt),g(nrmt),gp(nrmt,4)
C ... Local parameters
      integer konfig,nn,nre
      double precision dnu,eb1,eb2,pi,slo(5),sum,tol,val(5)
      pi = 4d0*datan(1d0)
      eb1 = -30
      eb2 = 20
      tol = 1d-10

      konfig = pnu
      nn = konfig-l-1
      dnu = dtan(pi*(0.5d0-pnu))
      val(1) = rmt
      slo(1) = dnu+1d0
      enu=-0.5d0

      call rseq(eb1,eb2,enu,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nrmt,
     .  nre)
      val(1) = val(1)/dsqrt(sum)
      slo(1) = slo(1)/dsqrt(sum)

C      call phidot(z,l,v,enu,a,b,rofi,nrmt,g,val,slo,tol,nn,gp,phi,dphi,
C     .  phip,dphip,p)

      call phidx(1,z,l,v,0d0,0d0,rofi,nrmt,2,tol,enu,val,slo,nn,g,gp,
     .  phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
C     dphip = (slo(2)-phip)/rmt


c|    write(stdo,200) l,enu,p,phi,dphi,phip,dphip
c|200 format(' PP',i2,'  e',f10.5,'  p',f10.5,'  bc',4f10.5)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine popta4 


      subroutine popta5(lmax,rtab,etab,itab,z,pl,rmax,rmt,nr,nrmt,
     .  rofi,psi,v,g,a,b,spid)

C- Write wave functions to plot file
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmax  :maximum l for a given site
Ci   rtab  :smoothing radii for wavefunction, each l
Ci   etab  :smoothed hankel energies for wavefunction, each l
Ci   itab  :1 if a wave function calculated, 0 if not
Ci   z     :nuclear charge
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   rmax  :muffin-tin radius, in a.u.
Ci   rmt   :muffin-tin radius
Ci   nr    :number of radial mesh points
Ci   nrmt  :number points to muffin-tin radius
Ci   rofi  :radial mesh points
Ci   psi   :wave function tabulated on the rofi mesh
Ci   v     :spherical potential (atomsr.f)
Ci   g     :normalized wave function times r
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   spid
Co Outputs
C    wave functions written to disk
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer itab(0:1),lmax,nr,nrmt,n0
      parameter (n0=10)
      double precision a,b,rmax,rmt,z,rtab(0:1),etab(0:1),
     .  rofi(nr),psi(nr,0:1),g(nr,2),v(nr),pl(0:n0-1)
      character spid*8
C ... Local parameters
      integer i,ifi,konfig,l,lgunit,lp1,m,n,nn,nre,fopna,stdo
      integer ltab(n0)
      double precision asm,dfl,drdi,eb1,eb2,eh,evl,fac,fl,flp1,r,rsm,
     .  slo,sum1,sum2,tol,val,wt,xi(0:20)
      character str*32

      stdo = lgunit(1)
      eb1 = -20
      eb2 = 20
      tol = 1d-8
      n = 0

      do  10  l = 0, lmax
        if (itab(l) .eq. 0) goto 10
        n = n+1
        ltab(n) = l
        lp1 = l+1
        rsm = rtab(l)
        eh = etab(l)
        asm = 1d0/rsm
        konfig = pl(l)
        nn = konfig-l-1

C ...   Smooth hankel fct outside rmt
        sum2 = 0d0
        do  12  i = nrmt, nr
          r = rofi(i)
          call hansmr(r,eh,asm,xi,l)
          psi(i,n) = xi(l)*(r**lp1)
          wt = 2*(mod(i+1,2)+1)/3d0
          if (i.eq.nrmt .or. i.eq.nr) wt=1d0/3d0
          drdi = a*(r+b)
          sum2 = sum2 + wt*drdi*psi(i,n)**2
   12   continue

C ...   Attach numerical solution inside MT sphere
        call hansmr(rmt,eh,asm,xi,l+1)
        fl   = xi(l)*rmt**l
        flp1 = xi(l+1)*rmt**(l+1)
        dfl  = l*fl/rmt-flp1
        val = rmt*fl
        slo = rmt*dfl+fl
        evl = -0.5d0
        call rseq(eb1,eb2,evl,tol,z,l,nn,val,slo,v,
     .     g,sum1,a,b,rofi,nrmt,nre)
        fac = val/(g(nrmt,1)*dsqrt(sum1+sum2))
        do  14  i = 1, nrmt
          psi(i,n) = fac*g(i,1)
   14   continue
        fac = 1d0/dsqrt(sum1+sum2)
        do  16  i = nrmt+1,nr
          psi(i,n) = psi(i,n)*fac
   16   continue
   10 continue

C ... Write the plot file
      write (str,'(''wfa_'',a)') spid
      write (stdo,344) str
  344 format(/' Write fit wavefunctions to plot file: ',a)
      ifi = fopna(str,-1,0)
      write (ifi,490) spid,rmax,rmt,(ltab(i),i=1,n)
  490 format('# Free-atom opt basis (divided by r) for species ',
     .   a/'# rmax=',f7.3,'   rmt=',f7.3,
     .   '   l=',8i2)
        write (ifi,'(''% rows '',i5,'' cols '',i3)') nr,n+1
      do  30  i = 1, nr
        write (ifi,495) rofi(i),(psi(i,m),m=1,n)
  495   format(f9.5,1p,8d14.5)
   30 continue
      call fclose(ifi)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine popta5 


      subroutine fctail(nr,nrmt,a,b,rsm,rofi,rhoc,c,eh)

C- Fit one Hankel to tail of core density.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nrmt  :number points to muffin-tin radius
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rsm   :smoothing radius
Ci   rofi  :radial mesh points
Ci   rhoc  :core density
Co Outputs
Co   c     :coefficient to fit of rhoc(spin+)+rhoc(spin-)
Co   eh    :energy
Cl Local variables
Cl   rmt   :muffin-tin radius
Cr Remarks
Cb Bugs
Cb   Should this be fit to smoothed function??
Cu Updates
Cu   19 Apr 02 Make rmt a local variable.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nrmt
      double precision a,b,c,eh,rmt,rsm,rofi(nr),rhoc(nr,2)
C ... Local parameters
      integer i,nsp,lgunit,stdo,ipr,nglob
      double precision ak1,akap,fit,q,q0,r,s,v0,wt
      character sout*80

      call getpr(ipr)
      stdo = lgunit(1)
      nsp =  nglob('nsp')
      rmt = rofi(nrmt)
      q0 = 0d0
      do  10  i = nrmt, nr
        r = rofi(i)
        wt = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
        if (i.eq.nrmt .or. i.eq.nr) wt = a*(r+b)/3d0
        q0 = q0 + wt*(rhoc(i,1)+rhoc(i,nsp))/(3-nsp)
   10 continue
      v0 = (rhoc(nrmt,1)+rhoc(nrmt,nsp))/(3-nsp)/(rmt*rmt)
      sout = ' '
      call awrit3('%?#(n>=30)#%N## coretail: q=%;3g, rho(rmt)=%;3g.',
     .  sout,len(sout),0,ipr,v0,q0)
C      if (ipr .ge. 20) write (stdo,339) v0,q0
C  339 format(/' coretail:  rho(rmt)=',f12.8,'   charge=',f12.8)
      if (dabs(q0) .lt. 1d-6) then
        c = 0d0
        eh = -1d0
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
      endif

C ... Parameters of hankel fct
      s = dsqrt(rmt**4 * v0**2 + 4*rmt*q0*v0)
      ak1 = (rmt*rmt*v0+s)/(2d0*q0)
C     ak2 = (rmt*rmt*v0-s)/(2d0*q0)
c|      write(stdo,975) ak1,ak2
c|  975 format('ak1,ak2=',2f14.8)
      akap = ak1
      c = rmt*v0*dexp(akap*rmt)
      eh = -akap*akap

      if (ipr .ge. 20) then
        call awrit2('%a  Fit with Hankel e=%;5g  coeff=%;5g',sout,
     .    len(sout),-stdo,eh,c)
      endif

C ... Test
      if (ipr .gt. 30) then
      write (stdo,501)
      q = 0d0
      do  20  i = nrmt, nr
        r = rofi(i)
        wt = 2*(mod(i+1,2)+1)*a*(r+b)/3d0
        if (i.eq.nrmt .or. i.eq.nr) wt = a*(r+b)/3d0
        fit = c*dexp(-akap*r)*r
        q = q+wt*fit
        if ((rhoc(i,1)+rhoc(i,nsp))/(3-nsp) .lt. 1d-8) goto 90
        if (mod(i,5).eq.0 .or. i.eq.nrmt)
     .     write (stdo,500) r,(rhoc(i,1)+rhoc(i,nsp))/(3-nsp),fit
  500   format(f12.6,2f14.8)
  501   format(6x,'r',12x,'rhoc',10x,'fit')
   20 continue
   90 continue
c|      v=c*dexp(-akap*rmt)/rmt
c|      write(stdo,885) q,q0,v,v0
c|  885 format('q,q0,v,v0=',4f14.8)
      endif

c ... look at smoothed core..
c|      rg=0.4
c|      qc=36
c|      sum0=-c*dexp(eh*rg*rg/4d0)/eh
c|      cg=qc-sum0
c|      write(stdo,888) qc,sum0,cg
c|  888 format(' qcore=',f10.4,'  sum0=',f12.6,'   cg=',f12.6)
c|      ag=1d0/rg
c|      fac=4d0*pi*(ag*ag/pi)**1.5d0
c|      q=0d0
c|      do i=1,nr
c|        r=rofi(i)
c|        wt=2*(mod(i+1,2)+1)*a*(r+b)/3d0
c|        if (i.eq.1 .or. i.eq.nr) wt=a*(r+b)/3d0
c|        call hansmr(r,eh,ag,xi,1)
c|        fit=c*xi(0)*r*r + cg*fac*dexp(-ag*ag*r*r)*r*r
c|        if (rhoc(i).gt.1d-10) write(49,490) r,rhoc(i),fit
c|  490   format(f12.6,2f16.8)
c|        q=q+wt*fit
c|      enddo
c|      write(stdo,965) q
c|  965 format(' integral over smoothed core:',f10.5)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine fctail 


