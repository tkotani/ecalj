Cgetarg...info...           structure ['asvsph2', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['asvsph2', 'array', 'sarray']
Cgetarg...info...           structure ['asvsph2', 'pot', 'spot']
Cgetarg...info...           structure ['asvsph2', 'ham', 'sham']
Cgetarg...info...           structure ['asvsph2', 'lat', 'slat']
Cgetarg...info...           structure ['asvsph2', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['asvsph2', 'm_struc_def']
      subroutine asvsph2(sctrl,slat,sarray,sspec,sham,spot,vne,rhos,
     .  nrhos,imake,ehterm,lhave)
       
       use m_struc_def  !Cgetarg

C- Make ASA sphere potential, potential parameters
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nbasp nclass nspec zbak nl nspin lncol loptc
Ci                lpgf lves lrel
Ci     Stored:
Ci     Passed to: lgors atscpp asamad
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: vol avw
Ci     Stored:
Ci     Passed to: asamad
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp ormax ohave oics oclabl onrcp oipc
Ci     Stored:
Ci     Passed to: asamad
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: eref
Ci     Stored:
Ci     Passed to: atscpp asamad
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: nmto kmto eterms
Ci     Stored:    eterms ehk thrpv seref amgm
Ci     Passed to:
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: orhrmx oves ovdif ovrmax ovintr opnu oqnu oqt opprel
Ci                opp osop ogrrme opmpol
Ci     Stored:    vmtz
Ci     Passed to: asamad
C     implicit none
C ... Passed parameters
      integer lhave,nrhos
Cgetarg       double precision sctrl(1),slat(1),sarray(1),sspec(1),sham(1),
Cgetarg      .  spot(1),ehterm,rhos(2,0:2,nrhos,2,2,1),vne
       real(8):: ehterm , rhos(2,0:2,nrhos,2,2,1) , vne 
       type(s_ctrl)::sctrl
       type(s_lat)::slat
       type(s_array)::sarray
       type(s_spec)::sspec(*)
       type(s_ham)::sham
       type(s_pot)::spot

c     ehterm(5)
C ... Local parameters
      logical bittst,lgors,lehk,havedc
      integer oclabl,ogrrme,oics,oinitc,oeula,oipc,onrcp,opmpol,opnu,
     .  opp,oqnu,oqt,orhrmx,ormax,osop,ovdif,oves,ovwk,ovintr,ovrmax,
     .  owk,opprel,oeulat,obxc
      integer ic,nclass,nclasp,nclspp,nbas,nl,nsp,igets,i,iclbsj,
     .  nbasp,lpgf,imake,is,lncol,loptc,lnsph,isw,lgunit,lves,k,ival,
     .  nrclas,nspec,neul
      integer mpipid,procid
      double precision thrpv,amgm,ehkk,emad,trumad,vmtz(2),zbak(2),
     .  dgetss,wk(5),dval,amag(3),dsqrt,ddot
C ... For atscpp
      character clabl*8
      double precision amgmat,amom,ehkat,sevat,thrpva,vol,avw,seref
C ... For third-generation LMTO
      integer nmto
      double precision kmto(20)
C     double precision qat,zat,dval
C ... Heap
      integer w(1)
      common /w/ w
C     This ordering must match sham->eterms; see uham
      double precision eterms(20),eterma(20)

      stop 'this routine should call mofr directly'

Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nbasp , integer(8) :: sctrl%nbasp 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg...info...           real(8) :: zbak(2) , real(8) :: sctrl%zbak 2
Cgetarg       call upack('ctrl nbas nbasp nclass nspec zbak',sctrl,
Cgetarg      .  nbas,nbasp,nclass,nspec,zbak)
       
       nbas=sctrl%nbas
       nbasp=sctrl%nbasp
       nclass=sctrl%nclass
       nspec=sctrl%nspec
       i_copy_size=size(sctrl%zbak) 
       call dcopy(i_copy_size,sctrl%zbak,1,zbak,1) 

Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg...info...           integer :: lncol , integer(8) :: sctrl%lncol 1
Cgetarg...info...           integer :: loptc , integer(8) :: sctrl%loptc 1
Cgetarg       call upack('ctrl nl nspin lncol loptc',sctrl,nl,nsp,lncol,loptc,0)
       
       nl=sctrl%nl
       nsp=sctrl%nspin
       lncol=sctrl%lncol
       loptc=sctrl%loptc

Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg...info...           integer :: ormax , integer(8) :: sarray%ormax 1
Cgetarg       call upack('array nclasp ormax',sarray,nclasp,ormax,0,0,0)
       
       nclasp=sarray%nclasp
       ormax=sarray%ormax

Cgetarg...info...           integer :: oinitc , integer(8) :: sarray%ohave 1
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: onrcp , integer(8) :: sarray%onrcp 1
Cgetarg...info...           integer :: oipc , integer(8) :: sarray%oipc 1
Cgetarg       call upack('array ohave oics oclabl onrcp oipc',sarray,oinitc,
Cgetarg      .  oics,oclabl,onrcp,oipc)
       
       oinitc=sarray%ohave
       oics=sarray%oics
       oclabl=sarray%oclabl
       onrcp=sarray%onrcp
       oipc=sarray%oipc

Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg...info...           integer :: ovdif , integer(8) :: spot%ovdif 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg...info...           integer :: ovintr , integer(8) :: spot%ovintr 1
Cgetarg       call upack('pot orhrmx oves ovdif ovrmax ovintr',spot,orhrmx,
Cgetarg      .  oves,ovdif,ovrmax,ovintr)
       
       orhrmx=spot%orhrmx
       oves=spot%oves
       ovdif=spot%ovdif
       ovrmax=spot%ovrmax
       ovintr=spot%ovintr

Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: oqt , integer(8) :: spot%oqt 1
Cgetarg...info...           integer :: opprel , integer(8) :: spot%opprel 1
Cgetarg       call upack('pot opnu oqnu oqt opprel',spot,opnu,oqnu,oqt,opprel,0)
       
       opnu=spot%opnu
       oqnu=spot%oqnu
       oqt=spot%oqt
       opprel=spot%opprel

Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
Cgetarg...info...           integer :: obxc , integer(8) :: spot%obxc 1
Cgetarg...info...           integer :: ogrrme , integer(8) :: spot%ogrrme 1
Cgetarg...info...           integer :: opmpol , integer(8) :: spot%opmpol 1
Cgetarg       call upack('pot opp osop obxc ogrrme opmpol',spot,
Cgetarg      .  opp,osop,obxc,ogrrme,opmpol)
       
       opp=spot%opp
       osop=spot%osop
       obxc=spot%obxc
       ogrrme=spot%ogrrme
       opmpol=spot%opmpol

Cgetarg...info...           integer :: nmto , integer(8) :: sham%nmto 1
Cgetarg...info...           real(8) :: kmto(20) , real(8) :: sham%kmto 6
Cgetarg       call upack('ham nmto kmto',sham,nmto,kmto,0,0,0)
       
       nmto=sham%nmto
       i_copy_size=size(sham%kmto) 
       call dcopy(i_copy_size,sham%kmto,1,kmto,1) 

Cgetarg...info...           integer :: oeula , integer(8) :: sham%oeula 1
Cgetarg...info...           integer :: neul , integer(8) :: sham%neula 1
Cgetarg       call upack('ham oeula neula',sham,oeula,neul,0,0,0)
       
       oeula=sham%oeula
       neul=sham%neula

Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg...info...           real(8) :: avw , real(8) :: slat%avw 1
Cgetarg       call upack('lat vol avw',slat,vol,avw,0,0,0)
       
       vol=slat%vol
       avw=slat%avw

Cgetarg       lpgf = igets('ctrl lpgf',sctrl)
       lpgf = int(sctrl%lpgf(1)) 

Cgetarg       lnsph = isw(lgors('ctrl lasa,32',sctrl))
       lnsph = isw ( iand(32,int(sctrl%lasa)) .ne.0 ) 

Cgetarg       lves = 2*igets('ctrl lves,1',sctrl)
       lves = 2 * iand(1,int(sctrl%lves)) 

      if (lpgf .ne. 0) then
        lves = lves+100
      endif

      ehkk = 0
      amgm = 0
      thrpv = 0
      seref = 0
      nclspp = 2*nclasp-nclass
C      call dpzero(ehterm,4)
C      call dpzero(amag,3)
Cgetarg...info...           real(8) :: eterms(20) , real(8) :: sham%eterms 20
Cgetarg       call upack1('ham eterms',sham,eterms) 
       
       i_copy_size=size(sham%eterms) 
       call dcopy(i_copy_size,sham%eterms,1,eterms,1) 

C     asvsph only modifies eterms(3:11)
      call dvset(eterms,3,11,0d0)
      call dvset(eterms,15,16,0d0)

C ... MPI: only master does sphere program
      procid = mpipid(1)
      if (procid .eq. 0) then

      call togprt()
      call defrr(owk,-nclspp)
      lehk = .true.
C      lhave = 1
      havedc = .true.
      do  40  ic = 1, nclasp
        nrclas = ival(w(onrcp),ic)
C   ... Classes to skip
        i = iclbsj(ic,w(oipc),-nbasp,1)
        if (i .lt. 1) goto 40
C   ... What to make
        if (lpgf .eq. 2 .and. i .le. nbas) goto 40
C       is = w(oics+ic-1)
        is = ival(w(oics),ic)

C       call pshpr(80)
        call dpzero(eterma,20)
Cgetarg...info...           integer :: oeula , integer(8) :: sham%oeula 1
Cgetarg...info...           integer :: neul , integer(8) :: sham%neula 1
Cgetarg         call upack('ham oeula neula',sham,oeula,neul,0,0,0)
         
         oeula=sham%oeula
         neul=sham%neula

        call defrr(oeula,nbas*neul*3)
        call atscpp2(sctrl,sspec,w(oeula),neul,nbas,is,ic,w(oclabl),
     .    w(ormax),imake,nl,nsp,w(oinitc+ic-1),zbak(1)/vol,avw,w(opnu),
     .    w(oqnu),w(oves),w(owk),0d0,w(opp),w(opprel),ehkat,sevat,
     .    w(oqt),amgmat,w(orhrmx),w(ovrmax),thrpva,w(osop),w(ogrrme),
     .    w(opmpol),w(ovintr),clabl,eterma,rhos,max(nl,neul))
        return
C       If estat energy changed, double-counting terms are available
        havedc = havedc .and. eterma(3) .ne. 0
        if (havedc) then
          eterma(4) = dval(w(oves),ic) * dval(w(oqt),ic) / 2
          call dpadd(eterms,eterma,3,11,dble(nrclas))
          call dpadd(eterms,eterma,15,16,dble(nrclas))
          call dpadd(eterms,eterma,19,19,dble(nrclas))
        else
          call dvset(eterms,3,11,0d0)
        endif
C       call poppr
Cgetarg         seref = seref + dgetss('spec eref',is,sspec)*nrclas
         seref = seref + (sspec(is)%eref) * nrclas 

        if (imake .eq. 4) goto 40

C   ... Check that potential and other pars available if to continue
        k = w(oinitc+ic-1)
Cgetarg         if (.not. lgors('ctrl lgen3',sctrl)) then
         if ( .not. iand(0,int(sctrl%lgen3)) .ne.0 ) then 

        if (mod(k/2,2) .ne. 1
     .    .or.(bittst(lncol,4).or.bittst(lncol,8)).and.(mod(k/4,2)).ne.1
     .    .or. loptc.gt.0 .and. (mod(k/32,2)) .ne. 1
     .    .or. lnsph .ne. 0 .and. (mod(k/16,2)) .ne. 1) then
          call awrit0(' LM (warning) class '//clabl//
     .     '%a missing pot, s-orbit or optical parms',' ',-80,lgunit(1))
C          lhave = 0
          goto 40
        endif
        endif

C   ... Only accumulate into ehkk if ehkat is set
C       qat = dval(w(oqt),ic)
C       zat = dgetss('spec z',is,sspec)
        if (havedc .and. lehk) then
Cgetarg           ehkat = ehkat - dgetss('spec eref',is,sspec)
           ehkat = ehkat - (sspec(is)%eref) 

          ehkk = ehkk + ehkat*nrclas
        else
          lehk = .false.
        endif
C        call atsev(ic,w(onrcp),w(oves),nl,nsp,w(oqnu),sevat,ehterm(2))
C        amgm = amgm + amgmat*nrclas
C        thrpv = thrpv + thrpva*nrclas
   40 continue
      call rlse(owk)
      call togprt()

C --- Make ves, emad, vmtz for new sphere Q, potential ---
C     imake=4: don't update pot->ves; use temporary array
      if (imake .eq. 4) then
        call defrr(ovwk,-nclspp)
      else
        ovwk = oves
      endif
      call asamad(sctrl,spot,sarray,slat,sspec,lves,
     .  w(opnu),w(oqnu),vne,w(ovwk),emad,trumad,vmtz)
C ... V[q] needed for total energy if a potential shift
      call dpzero(w(ovdif),nclspp)
      if (mod(lves,100) .eq. 2) then
        call pshpr(0)
        call asamad(sctrl,spot,sarray,slat,sspec,0,w(opnu),w(oqnu),
     .    vne,w(ovdif),wk(1),wk(2),wk(3))
        call daxpy(nclspp,-1d0,w(ovwk),1,w(ovdif),1)
        call poppr
      endif
      if (imake .eq. 4) call rlse(ovwk)

C ... Add Madelung contribution to total electrostatic energy
      eterms(3) = eterms(3) + emad
C ... Total magnetization, noncollinear case
      if (neul .gt. 0) then
        amgm = dsqrt(ddot(3,amag,1,amag,1))
        eterms(15) = amgm
      endif
C ... Repack double-counting terms
Cgetarg...info...           real(8) :: eterms(20) , real(8) :: sham%eterms 20
Cgetarg       call pack1('ham eterms',sham,eterms)
       
       i_copy_size=size(sham%eterms) 
       call dcopy(i_copy_size,eterms,1,sham%eterms,1) 


C --- Make NMTO potential parameters ---
C      if (lgors('ctrl lgen3',sctrl)) then
C        call ppmesh(sctrl,sspec,slat,nl,nmto,kmto,w(oclabl),
C     .    vmtz,w(oics),w(oppn))
C      endif

C --- Save energy for HF (to recover if information lost) ---
C terms:
C (1) ehkk (no Mad) (2) sumeV(VH=0)
C (3) sum Q_i V(R) (4) emad
C (5) d.c. terms from applied fields.
C      ehterm(4) = emad
C      if (lehk) then
C        ehterm(1) = ehkk
C        ehkk = ehkk + emad - ehterm(5)
C      endif

      thrpv = thrpv + trumad
Cgetarg...info...           real(8) :: vmtz(2) , real(8) :: spot%vmtz 1
Cgetarg       call pack1('pot vmtz',spot,vmtz)
       
       spot%vmtz=vmtz(1)

Cgetarg...info...           real(8) :: ehkk , real(8) :: sham%ehk 1
Cgetarg...info...           real(8) :: thrpv , real(8) :: sham%thrpv 1
Cgetarg...info...           real(8) :: seref , real(8) :: sham%seref 1
Cgetarg...info...           real(8) :: amgm , real(8) :: sham%amgm 1
Cgetarg       call pack5('ham ehk thrpv seref amgm',sham,ehkk,thrpv,seref,amgm,
Cgetarg      .  0)
       
       sham%ehk=ehkk 
       sham%thrpv=thrpv 
       sham%seref=seref 
       sham%amgm=amgm 

C      if (lehk) lhave = lhave+10

      endif

C --- MPI broadcast everything passed out of asvsph ---
C     sarray and contents of sarray pointers
      call mpibc1(w(oinitc),nclasp,2,.false.,'asvsph','initc')
C     sham
      call mpibc1(sham,nint(sham%size)*nbas,4,.false.,'asvsph','sham')
C     spot and and contents of spot pointers
      call mpibc1(w(orhrmx),nclasp,4,.false.,'asvsph','rhrmx')
      call mpibc1(w(oves),nclspp,4,.false.,'asvsph','ves')
      call mpibc1(w(ovdif),nclspp,4,.false.,'asvsph','vdif')
      call mpibc1(w(ovrmax),2*nclspp,4,.false.,'asvsph','vrmax')
      if (ovintr .gt. 0) then
        i = nclasp*(nl*nsp)**2
        call mpibc1(w(ovintr),i,4,.false.,'asvsph','vintr')
      endif
      call mpibc1(w(opnu),nl*nsp*nclasp,4,.false.,'asvsph','pnu')
      call mpibc1(w(oqnu),3*nl*nsp*nclasp,4,.false.,'asvsph','qnu')
      call mpibc1(w(oqt),nclasp,4,.false.,'asvsph','qt')
      call mpibc1(w(opp),6*nl*nsp*nclasp,4,.false.,'asvsph','pp')
C     lrel = igets('ctrl lrel',sctrl)
C     if (lrel .eq. 2) then
      if (opprel .gt. 0) then
        call mpibc1(w(opprel),32*nl*nl*max(nclspp,nspec),4,.false.,
     .    'asvsph','pprel')
      endif
C     if (igets('ctrl lncol',sctrl) .gt. 0) then
      if (osop .gt. 0) then
       call mpibc1(w(osop),9*nl*nsp*nsp*nclasp,4,.false.,'asvsph','sop')
      endif
C     if (igets('ctrl loptc',sctrl) .gt. 0) then
      if (ogrrme .gt. 0) then
       call mpibc1(w(ogrrme),8*nl*nsp*nclasp,4,.false.,'asvsph','grrme')
      endif
C     if (lgors('ctrl lasa,32',sctrl)) then
      if (opmpol .gt. 0) then
       i = (2*nl-1)*nl**2*3*nsp*nclasp
       call mpibc1(w(opmpol),i,4,.false.,'asvsph','pmpol')
      endif
      call mpibc1(spot,nint(spot%size),4,.false.,'asvsph','spot')
C     Remaining arguments passed to asvsph
C      call mpibc1(ehterm,5,4,.false.,'asvsph','ehterm')
C      call mpibc1(lhave,1,2,.false.,'asvsph','lhave')

Cgetarg       end
       end subroutine asvsph2 


      subroutine mofr(x1,x2,x3,mag,sdir,yl,rsq,sqmag,iofr,z,lmax,nl,a,
     .  b,nr,rofi,v,tol,qnu,eula,neula,nbas,rmax,nsp,pnu,npts,nedge,
     .  rhos,nrhos)

C-
C ----------------------------------------------------------------------
Ci Inputs
Ci   x1
Ci   x2
Ci   x3
Ci   mag
Ci   sdir
Ci   yl
Ci   rsq
Ci   sqmag
Ci   iofr
Ci   z     :nuclear charge
Ci   lmax  :maximum l for a given site
Ci   nl    :(global maximum l) + 1
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   tol
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   eula  :Euler angles for noncollinear spins
Ci   neula
Ci   nbas  :size of basis
Ci   rmax  :augmentation radius, in a.u.,
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   npts  :number of DOS tabulation points (input; for sampling only)
Ci   nedge
Ci   rhos  :spin density-matrix
Ci   nrhos :number of channels which spin density-matrix is stored
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   16 Sep 04 (S. Brink) first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,neula,nr,nl,npts,lmax,iofr(npts),nrhos,nedge
      double precision a,b,eula(nbas,neula,3),mag(nbas,npts,3),
     .  pnu(nl,nsp),qnu(3,nl,nsp),rhos(2,0:2,nrhos,2,2,1),rmax,
     .  rofi(nr,2),rsq(npts),sdir(nbas,neula,3),sqmag(npts),tol,
     .  v(nr,nsp),x1(npts),x2(npts),x3(npts),yl(npts,nl**2),z
C ... Local parameters
      integer fopna,i,j,k,l,ifi,isp,lm,low,m,nn,nre,intrvl,pt
      double precision eb1,eb2,eval,fac,g(nr,2),pi,slo,sum,val,wavfn
C     double precision sumrho,summag(3)
C ... Heap
      integer w(1)
      common /w/ w

C     stdo = lgunit(1)
      pi = 4d0*datan(1d0)
C     sumrho = 0d0
C     summag(1) = 0d0
C     summag(2) = 0d0
C     summag(3) = 0d0
      intrvl = 2*rmax/(nedge-1)
      do  i = 1, npts
        j = (i-1)/nedge**2
        k = (i-1-j*nedge**2)/nedge
        x1(i) = -rmax + j*intrvl
        x2(i) = -rmax + k*intrvl
        x3(i) = -rmax + (i-1-j*nedge**2-k*nedge)*intrvl
      end do
      do  i = 1, nbas
        do  lm = 1, neula
          sdir(i,lm,1) = dsin(eula(i,lm,2))*dcos(eula(i,lm,1))
          sdir(i,lm,2) = dsin(eula(i,lm,2))*dsin(eula(i,lm,1))
          sdir(i,lm,3) = dcos(eula(i,lm,2))
        end do
      end do
      do  i = 1, npts
        sqmag(i) = x1(i)**2+x2(i)**2+x3(i)**2
        if ( sqmag(i) .eq. 0.0d0 ) then
          sqmag(i) = 1
        end if
        if ( sqmag(i) .gt. rmax**2 ) then
          sqmag(i) = 100
        else
          x1(i) = x1(i)/(dsqrt(sqmag(i)))
          x2(i) = x2(i)/(dsqrt(sqmag(i)))
          x3(i) = x3(i)/(dsqrt(sqmag(i)))
        end if
      end do
      call ropyln(npts,x1,x2,x3,lmax,npts,yl,rsq)
      call dpzero(mag,nbas*npts*3)
      do  i = 1, nbas
        do  isp = 1, nsp
          fac = 2*(1.5d0-dble(isp))
          do  l = 0, lmax
            nn = int(pnu(l+1,isp)) - l - 1
            val = rmax
            slo = 1 + dtan(pi*(0.5d0 - pnu(l+1,isp)))
            eb1 = -50d0
            eb2 =  20d0
            eval= -15d0
            call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v(1,isp),g,
     .          sum,a,b,rofi,nr,nre)
            do  pt = 1, npts
              if (sqmag(pt) .ne. 100d0) then
              do m = -l, l
                lm = l**2+l+m+1
                low = aint((dlog(dsqrt(sqmag(pt))/b+1))/a)+1
                wavfn=(g(low,1)+(g(low+1,1)-g(low,1))*(dsqrt(sqmag(pt))-
     .          rofi(low,1))/(rofi(low+1,1)-rofi(low,1)))*yl(pt,lm)
                if ( neula .eq. nl**2 ) then
                  mag(i,pt,1) = mag(i,pt,1) + wavfn**2*fac*
     .            rhos(1,0,lm,isp,isp,1)*sdir(i,lm,1)
                  mag(i,pt,2) = mag(i,pt,2) + wavfn**2*fac*
     .            rhos(1,0,lm,isp,isp,1)*sdir(i,lm,2)
                  mag(i,pt,3) = mag(i,pt,3) + wavfn**2*fac*
     .            rhos(1,0,lm,isp,isp,1)*sdir(i,lm,3)
                elseif ( neula .eq. nl ) then
                  mag(i,pt,1) = mag(i,pt,1) + wavfn**2*fac*
     .            rhos(1,0,l+1,isp,isp,1)*sdir(i,l+1,1)/(2*l+1)
                  mag(i,pt,2) = mag(i,pt,2) + wavfn**2*fac*
     .            rhos(1,0,l+1,isp,isp,1)*sdir(i,l+1,2)/(2*l+1)
                  mag(i,pt,3) = mag(i,pt,3) + wavfn**2*fac*
     .            rhos(1,0,l+1,isp,isp,1)*sdir(i,l+1,3)/(2*l+1)
                else
                  mag(i,pt,1) = mag(i,pt,1) + wavfn**2*fac*
     .            rhos(1,0,l+1,isp,isp,1)*sdir(i,1,1)/(2*l+1)
                  mag(i,pt,2) = mag(i,pt,2) + wavfn**2*fac*
     .            rhos(1,0,l+1,isp,isp,1)*sdir(i,1,2)/(2*l+1)
                  mag(i,pt,3) = mag(i,pt,3) + wavfn**2*fac*
     .            rhos(1,0,l+1,isp,isp,1)*sdir(i,1,3)/(2*l+1)
                end if
              end do
            else
              mag(i,pt,1) = 0d0
              mag(i,pt,2) = 0d0
              mag(i,pt,3) = 0d0
            end if
            end do
          end do
        end do
      end do
      ifi = fopna('mofr',-1,0)
      write(ifi,'(2x,"x",10x,"y",10x,"z",10x,"mx",10x,"my",10x,"mz")')
      do  i = 1, nbas
        do  pt = 1, npts
          if (sqmag(pt) .ne. 100d0) then
            x1(pt) = x1(pt)*(dsqrt(sqmag(pt)))
            x2(pt) = x2(pt)*(dsqrt(sqmag(pt)))
            x3(pt) = x3(pt)*(dsqrt(sqmag(pt)))
          end if
c           if ((x1(pt).ne.0d0).and.(x2(pt).ne.0d0).and.(x3(pt).ne.0d0))
c     .       then
          write(ifi,'(f8.4,3x,f8.4,3x,f8.4,3x,f8.4,3x,f8.4,3x,f8.4)')
     .      x1(pt),x2(pt),x3(pt),mag(i,pt,1),mag(i,pt,2),mag(i,pt,3)
c          endif
        end do
      end do
      call fclose(ifi)
Cgetarg       end
       end subroutine mofr 


