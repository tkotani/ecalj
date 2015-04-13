#define USEOPTIONAL 1
!! ------------------------------
!! ------------------------------
      subroutine sugw ( ssite , sspec , slat , sham , nbas , ndham 
     .      , smpot , vconst , lcplxp , sv_p_osig , sv_p_otau , sv_p_oppi 
     .      , ppn , vrmt , spotx , sv_p_osigx , sv_p_otaux , sv_p_oppix , 
     .      jobgw ,nk1,nk2,nk3,sigmamode,ndimsig, numproc,procid,master)
      use m_struc_def,only:  s_rv1, s_lat, s_ham, s_spec, s_site
      use m_globalvariables,only: globalvariables
      use m_phmbls,only: phmbls
      implicit none
!! == Driver for fpgw (to prepare eigenfuncitons for fpgw) ==
!! NOTE: following documents are not carefully examined. Not believe everything.
!!i Inputs
!!i   ssite,sspec,slat,sham :struct defined in m_struc_def.
!!i   osig,otau,oppi  augmentation matrices, s_rv1
!!i   sham%rv_a_ohrs: real space Sigma_vxc 
!!i   nbas  :size of basis
!!i   smpot :smooth potential on uniform mesh (mkpot.f)
!!i   vconst:constant to be added to potential
!!i   lcplxp:0 if ppi is real; 1 if ppi is complex
!!i   ppn   :potential parameters, nmto style
!!i   vrmt  :electrostatic potential at MT boundaries
!!i   jobgw :-999 prompt for and read jobgw from stdin
!!i         :0 create files SYMOPS,LATTC,CLASS,NLAindx
!!i         :1 create files gwb,gw1,gw2,gwa,vxc,evec,rhoMT.*,normchk
!!o Outputs
!!o   Files written according to jobgw
!!o   The following shows what files are written for jobgw=1
!!o   and the records stored in each file.
!!o   gw1:  evals of Hlda+sigma-vxc
!!o        *for each q-point and spin:
!!o         q, evl(i:ndimh)
!!o   gw2:  evals of Hlda+sigma-vxc-vxc
!!o        *for each q-point and spin
!!o         q, evl(i:ndimh)
!!o   gwb:  Information about eigenfunctions, matrix elements
!!o         nat,nsp,ndima,ndham,alat,qlat,ef0,nqbz,plat,nqnum
!!o         lmxa(1:nat), bas(1:3,1:nat)
!!o         ngpmx,ngcmx  -- largest number of G vectors for psi, vcoul
!!o        *for each q-point and spin:
!!o           q, ndimh
!!o           evec, evl, cphi, ngp, ngc
!!o           ngvecp, ngvecc, pwz  (G vectors for psi,vcou; PW expansion of z)
!!o   gwa:  site data.
!!o        *for each site:
!!o           z, nr, a, b, rmt, lmaxa, nsp, ncore
!!o           konfig(0:lmaxa) : note
!!o           rofi: radial mesh
!!o          *for each l, spin isp
!!o             l, isp
!!o             radial w.f. gval: phi
!!o             radial w.f. gval: phidot
!!o             radial w.f. gval: phiz    written if konfig(l)>10
!!o             *for each l, spin, konf
!!o                icore, l, isp, konf, ecore(icore)+vshft
!!o                gcore(1:nr,1,icore)
!!o   evec: eigenvectors.
!!o         ndham, nsp, nnn, nqnum
!!o        *for each q-point and spin:
!!o           q, evec(1:ndimh,1:ndimh)
!!o   vxc:  matrix elements of XC potential
!!o         ndham, nsp, nnn
!!o        *for each q-point and spin:
!!o           q, vxc
!! --------------------------------------------------------------------------
Cl Local variables
Cl   lsig  :switch to create files vxc and evec for making sigma
Cl   lwvxc :T, write to evec and vxc files
Cl   nnn   :number of qp in IBZ
Cl   nqnum :number of qp in IBZ plus points related to Gamma. See mkqg.F to generate QGpsi
Cl         :or generally the number of qp at which eigenfunctions calc.
Cl   ngp   :no. G vectors for eigenfunction expansion (depends on q-pt)
Cl   ngc   :no. G vectors for coulomb interaction (depends on q-pt)
Cl   ispc  :2 when working on (2,2) block of noncollinear hamiltonian;
Cl         :otherwise 1
Cl   ipb   :index to true basis (excluding floating orbitals)
Cl         :given site index including those orbitals
Cl   ndima :number of augmentation channels
Cl   loldpw:0 use version 6 convention in calling pwmat
Cl         :1 call pwmat2 instead of pwmat
Cl         :2 call pwmat (or pwmat2) with shortened q vector
Cl         :3 both options 1 and 2
Cl         :See Remarks
Cb Bugs
Cb   code writes extra file evec whose data should be
Cb   extracted from gwb.
Cr Remarks
Cr   In version 6, the sugw passed an unshortened q to pwmat.
Cr   This cannot be used with an APW basis.
Cr   The 2's bit in loldpw controls which q is passed to pwmat.
Cr   The call to pwmat can optionally be replaced by a call to the
Cr   simpler and cleaner pwmat2 (see comments preceding the call to
Cr   pwmat or pwmat2).  However, the results are not identical.
Cr   The pwmat construction depends on both the LMTO cutoff gmax
Cr   and the GW cutoff QpGcut_psi, even though these nominally serve
Cr   the same purpose.
Cr   The pwmat2 construction depends only on the GW cutoff QpGcut_psi.
Cr   The pwmat construction is usually more accurate (i.e. the norm of
Cr   the overlap matrix is closer to 1; see output in file normchk),
Cr   because gmax is typically larger than QpGcut_psi.
Cu Updates
Cu   29 Jan 09 Incorporate APW basis
Cu   27 Mar 07 bug fix: expunging floating orbitals from class list, file CLASS
Cu   30 Aug 05 sugw handles ngp=0 and/or ngc=0
Cu    5 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu             Bug fix, job 5 case
Cu    4 Sep 04 Adapted to extended local orbitals
Cu    1 Sep 04 Adapted to handle complex ppi; S.O. put into ppi
Cu      Mar 04 (mark) small changes for bandmode
Cu      Sep 03 (takao)
Cu             Implemented job 5
Cu             read nqbz from QGpsi.  GWIN0 no longer used.
Cu   14 Sep 03 dimensioning bug fix when inequivalent lmxa
Cu   07 Jun 03 sugw redesigned for new interpolation mode for sigma
Cu             jobs 1 and 2 have been combined.
Cu             Altered argument list.
Cu   20 Jun 02 (S. Faleev) write vxc to disk
Cu   18 Jun 02 Added debugging code to check call to roth
Cu   25 Apr 02 Added local orbitals
Cu   25 Oct 01 (T. Kotani) split functions into parts, with job
Cu             Binary files imcompatible with prior versions.
Cu   23 Apr 01 First created
!! ----------------------------------------------------------------------
      integer:: i_copy_size,i_spackv
      integer nbas,ndham,n0,nkap0,nppn,lcplxp
      parameter (n0=10, nppn=12, nkap0=3)
      integer:: jobgw , lh(n0)
      type(s_rv1) :: sv_p_oppi(3,nbas)
      type(s_rv1) :: sv_p_otau(3,nbas)
      type(s_rv1) :: sv_p_osig(3,nbas)
      type(s_rv1) :: sv_p_oppix(3,nbas)
      type(s_rv1) :: sv_p_otaux(3,nbas)
      type(s_rv1) :: sv_p_osigx(3,nbas)
      type(s_site):: ssite(nbas)
      type(s_spec):: sspec(*)
      type(s_lat):: slat
      type(s_ham):: sham
      real(8):: rsml(n0), ehl(n0), ppn(nppn,n0,nbas), vconst, smpot(*), spotx(*), vrmt(nbas)
      logical :: bandmode=.false.,endofline=.false.,lwvxc
      logical cmdopt
      integer fopna,fopnx,i,i1,i2,iat,ib,ibr,icore,ierr,ifeigen,!ifi,
     .ifiqg,iflband(2),ifqeigen,ifsyml,igets,igetss,iix,iline, 
     .im1,im2,ipb(nbas),ipqn,ipr,iprint,iq,is,isp,ispc,j,job,jobb,k1,
     .k2,k3,konf,konfig(0:n0),l,lchk,ldim,loldpw, ldham(8,2),
     .lgunit,lmaxa,lmxax,lpdiag,lrsig,lsig,mx,mxint,n1,n2,n3,nat,
     .ncore,ndima,ndimh,nevl,nev,nglob,ngp,ngp_p, 
     .ngpmx,nline,nlinemax,nlmax,nmx,nn1,nn2,nnn,npgrp,
     .nphimx,npqn,nqbz,nqibz,nqnum,nqnumx,nqtot,nr,nsgrp,nsp,nspc,iqibz,imx, stdo,
     & ifigwb,ifigwa,ifinormchk,ifinlaindx,ifigw1,ifildima,ifigwn,
     & ificlass,ifievec,ifievecx,ifigw2,ifiqbz,ifievv
      complex(8) ,allocatable :: aus_zv(:)
      real(8) ,allocatable :: ww_rv(:)
      real(8) q(3),QpGcut_psi,QpGcut_cou,dum,dval,
     .xx(5),gmax,pnu(n0,2),pnz(n0,2),ecore(50),a,z,rmt(nbas),b,vshft,
     .alat,alfa,ef0,plat(3,3),qlat(3,3),qp(3),qpos,q_p(3), !,qx(3)
     .epsovl,dgets
      integer ,allocatable:: ips(:),ipc(:),ipcx(:),lmxa(:),
     .nlindx(:,:,:),ngvecp(:,:),ngvecp_p(:,:) !takao feb2012 ,ngvecc(:,:)
      integer,allocatable :: konft(:,:,:),iiyf(:),ibidx(:,:),nqq(:)
      real(8) ,allocatable:: wk(:,:),
     .bas(:,:),rofi(:),rwgt(:),gcore(:,:,:),gval(:,:,:,:,:),evl(:,:)
      real(8),allocatable:: ovv(:,:),evl_p(:,:),
     .qq1(:,:),qq2(:,:),cphin(:,:,:)
      complex(8),allocatable:: ham(:,:),ovl(:,:),evec(:,:),vxc(:,:),
     .ppovl(:,:),phovl(:,:),pwh(:,:),pwz(:,:),pzovl(:,:),
     .testc(:,:),testcd(:),ppovld(:),cphi(:,:,:),cphi_p(:,:,:),
     .geig(:,:,:),geig_p(:,:,:),sene(:,:)
      integer::isize_ham(3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
      integer ndhamx,nspx,nnnx,ifiv
      character strn*120
#if SIGNN
      complex(8),allocatable:: sig(:,:)
#endif
C ... For PW basis
c      integer:: oigv2
      integer ,allocatable :: kv_iv(:)
      integer pwmode,napw
      double precision pwemin,pwemax,pwgmin,pwgmax,eomin

C ... for reading self-energy
      integer nqsig
ckino Dec.14.2011: real(8),pointer :: rv_p_oqsig(:) =>NULL()

C ... for band plotting
      real(8),allocatable :: ovvx(:,:,:)
      real(8) ::ovvpp
      integer idxdn(n0,nkap0)
      character lsym(0:n0-1)*1, lorb(3)*1, dig(9)*1, strn4*4
      data lsym /'s','p','d','f','g','5','6','7','8','9'/
      data lorb /'p','d','l'/
      data dig /'1','2','3','4','5','6','7','8','9'/
C      integer jfi
C      integer,allocatable::  ngvecp_(:,:),ngvecc_(:,:)
C      real(8),allocatable:: evl_(:,:)
C      complex(8),allocatable:: evec_(:,:),cphi_(:,:,:),pwz_(:,:)

      integer:: w(1) !! dummy
      real(8):: dnn(3),qlatinv(3,3),qout(3),qtarget(3),qrr(3),axx,bxx,qxxx(3),qxx1(3),qxx2(3)
      integer:: inn(3),iqzz,nqzz,iiiii
      logical:: debug=.false.
      complex(8),allocatable:: evecout(:,:),evecr(:,:)
      real(8),allocatable:: qzz(:,:)
      integer,allocatable:: igv2(:,:)
      logical:: l_dummy_isanrg,isanrg,oncewrite,newsigmasw,sigmamode !,noshorbz

      integer(4):: nk1,nk2,nk3,   irr,nqi, nmcore !feb2012takao
c      complex(8):: sfz(nk1,nk2,nk3,ndham,ndham,*) !takao assume ndham =ndimh, OK?

      integer:: ndimsig
      
      integer:: procid,master,nbalance,numproc,iqq,nrr,ifiproc,iproc,iadd,ldw
      logical:: rank0
      integer,allocatable:: iprocq(:,:)
      character*256:: ext

      logical ::nexist,magexist  !! june2013 magfield is added for fsmom mode.
      real(8):: vnow
      integer:: ifimag,ifile_handle !=9078

      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx),iapw,ival
      real(8)::ppin(3)
      real(8):: rlatp(3,3),xmx2(3)
      real(8):: qqq(3)
      integer:: hreal,mode,iwdummy,lso,jx
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),senex(:,:),ovlmtoi(:,:),ovliovl(:,:) 

#if MPI | MPIK
      include "mpif.h"
#endif

!! --- Setup ---
      call tcn ('sugw')
      if(iprint()>10) print *,'sugw:'

!! june2013 magfield is added
      ifimag = ifile_handle()
      inquire(unit=ifimag,opened=nexist)
      if(nexist) stop 'sugw:ifimag used already (too sloppy manner ; fix this). '
      open(ifimag,file='MagField',status='old',err=112)
      magexist=.true.
      if(magexist) then
         open(ifimag,file='MagField',status='old')
         read(ifimag,*)vnow
         close(ifimag)
         write(6,"('Add mag.field to eval. -vnow/2 for isp=1, +vnow/2 for isp=2. vnow=',d13.6)")vnow
         goto 113
      endif   
 112  continue
      magexist=.false.
 113  continue

      rank0=.false.
      if(procid==master) rank0=.true.
      call getpr(ipr)
      alat=slat%alat
      n1 = slat%nabc(1)
      n2 = slat%nabc(2)
      n3 = slat%nabc(3)
      plat =slat%plat
      qlat =slat%qlat
      gmax =slat%gmax
      npgrp =slat%npgrp
      nsgrp =slat%nsgrp
      call fftz30(n1,n2,n3,k1,k2,k3)
      stdo = lgunit(1)
      lchk = 1
      nsp  = globalvariables%nsp
      nspc = globalvariables%nspc
      lrsig = int(sham%lsig)
      lwvxc = .not. cmdopt('--novxc',7,0,strn)
      nqsig = sham%nqsig

C     for now
      nphimx = 3

!!  Count number of atoms : exclude floating orbitals
      nat = 0
      do  i = 1, nbas
        is = ssite(i)%spec
        lmaxa = int(sspec(is)%lmxa)
        if (lmaxa .gt. -1) then
          nat = nat + 1
        endif
        ipb(i) = nat
      enddo

!! Select job
      job = jobgw
      if (job .eq. -999) then
        write(stdo,*) ' lmfgw: input one of the following jobs:'
        write(stdo,*) '  -1 : creates files:  GWinput, QPNT, QIBZ, Q0P, QGpsi, QGcou, KPTin1BZ'
        write(stdo,*) '   0 : init mode; creates files: SYMOPS, LATTC, CLASS, NLAindx, ldima'
        write(stdo,*) '   1 : GW setup mode; creates: files gwb,gw1,gw2,gwa,vxc,evec,rhoMT.*,normchk'
        write(stdo,*) '   4 : band mode '
        write(stdo,*) '   5 : eigenvalue-only mode ' !takao Sep 2003
        write(stdo,*) ' job?'
        read (5,*) job
      endif
      if(iprint()>10) print '(/a,i2)', ' gw setup, job', job
      l_dummy_isanrg=isanrg(job,0,5,'sugw:','job',.true.)
      if (job .eq. 3) call rxi('sugw: bad job = ',job)
      bandmode = .false.
      jobb = 0
      if (job .eq. 4) then
        job = 1
        bandmode = .true.
      elseif (job .eq. 5) then
        job = 1
        bandmode = .true.
        jobb = 5
      endif
      lsig = 1
      if (bandmode) lsig=0
      napw = 0
      i_copy_size=size(sham%ldham)
      call i8icopy(i_copy_size,sham%ldham,1,ldham,1)
      ldim=ldham(1,1)
      pwmode= sham%pwmode
      pwemin= sham%pwemin
      pwemax= sham%pwemax
      allocate(evl(ndham,nsp))
!! ndima
      ndima = 0
      lmxax = -1
      do  ib = 1, nbas
        is=ssite(ib)%spec
        lmaxa=sspec(is)%lmxa
        i_copy_size=size(sspec(is)%pz)
        call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1)
        lmxax = max(lmxax,lmaxa)
        if (lmaxa .gt. -1) then
          do  l = 0, lmaxa
            npqn = 2
            if (pnz(l+1,1) .ne. 0) npqn = 3
            ndima = ndima + npqn*(2*l+1)
          enddo
        endif
      enddo

!! === Make files SYMOPS, LATTC, CLASS, NLAindx, ldima ===
      if (job==0 .and. rank0) then
        call info(30,1,1,' Creating files SYMOPS, LATTC, CLASS, NLAindx, ldima',0,0)
        allocate(lmxa(nbas),bas(3,nbas)) 
!! ldima
        ifildima = fopnx('ldima',2,2,-1)
        rewind ifildima
        do  ib = 1, nbas
          lmxa(ib) = sspec(int(ssite(ib)%spec))%lmxa  !int(sspec(ips(ib))%lmxa)
          if (lmxa(ib) .gt. -1) then
            call orbl ( ib , 0 , ldim , sham%iv_a_oindxo , norb , ltab , ktab , xx , offl , i1 )
            write(ifildima,"(3i10)") i1
          endif
        enddo
!! LATTC
        lmxax = mxint(nbas,lmxa)
        allocate(konft(0:lmxax,nbas,nsp))
        do ib = 1, nbas
           is = ssite(ib)%spec
           call dcopy(size(ssite(ib)%pnu), ssite(ib)%pnu,1,pnu,1)
           call dcopy(size(ssite(ib)%pz),  ssite(ib)%pz, 1,pnz,1)
           do  isp = 1, nsp
              do  l  = 0, lmxa(ib)
                 konft(l,ib,isp) = pnu(l+1,isp)
                 if( mod(pnz(l+1,isp),10d0) .lt. pnu(l+1,isp) .and. pnz(l+1,isp) .gt. 0) then
                    konft(l,ib,isp) = mod(pnz(l+1,isp),10d0)
                 endif    
              enddo
           enddo
        enddo
        call wlattc(alat,plat,nbas,nat,ipb,lmxax,lmxa,nsp,konft)
!! NLAindx
        ifinlaindx = fopnx('NLAindx',2,2,-1)
        rewind ifinlaindx
        write(ifinlaindx,'(''----NLAindx start---------------''/I6)') ndima
        ndima = 0
!!    This loop order is backwardly compatible with prior versions
        do  ipqn = 1, 3
          do  ib = 1, nbas
            is =  ssite(ib)%spec
            lmaxa=sspec(is)%lmxa
            i_copy_size=size(ssite(ib)%pnu)
            call dcopy  (i_copy_size ,ssite(ib)%pnu,1,pnu,1)
            i_copy_size=size(ssite(ib)%pnu)
            call dcopy  (i_copy_size,ssite(ib)%pz, 1,pnz,1)
            call i8icopy(sspec(is)%idxdn,sspec(is)%idxdn,1,idxdn,1)
            if (lmaxa .gt. -1) then
              do  l = 0, lmaxa
                npqn = 2
                if (pnz(l+1,1) .ne. 0) npqn = 3
                if (ipqn .le. npqn) then
                  konf = pnu(l+1,1)
                  if (ipqn .eq. 3) konf = mod(pnz(l+1,1),10d0)
                  strn4 = dig(konf)//lsym(l)//'_'//lorb(ipqn)
                  if (idxdn(l+1,1) .eq. 1 .or. idxdn(l+1,2) .eq. 1)
     .            call chcase(0,1,strn4(2:2))
                  write(ifinlaindx,'(i6,i3,i4,i6,4x,a)')
     .            ipqn, l, ipb(ib), ndima, strn4
                  ndima = ndima + (2*l+1)
                endif
              enddo
            endif
          enddo
        enddo
        call fclr(' ',ifinlaindx)
!! CLASS
        ificlass = fopnx('CLASS',2,2,-1)
        rewind ificlass
        allocate(ipc(nbas),ipcx(nbas))
        do ib=1,nbas
           ipc(ib) = ssite(ib)%class !call spackv_array_copy_i8_i('u',ssite(ib)%class,1,ib,ipc)
        enddo
        call pvsug1(nbas,lmxa,ipc,ipcx)
        do  i = 1, nbas
           if(lmxa(i) .gt. -1) then
              write(ificlass,'(2I4)') ipb(i), ipcx(i)
           endif
        enddo
        deallocate(ipc,ipcx)
        call fclr(' ',ificlass)
        call fexit(0,1,' OK! '//'lmfgw mode=0 generated LATTC CLASS NLAindx ldima',0)
      endif

!! ==== Read file NLAindx ====
      allocate(nlindx(3,0:lmxax,nat))
      ifinlaindx = fopnx('NLAindx',2,1,-1)
      call ioaindx(3,lmxax,nat,ndima,nlindx,ifinlaindx)
!! ... Read QGpsi !takao feb2012 ---> Not QGcou now
      if (jobb .eq. 0) then
        ifiqg  = fopnx('QGpsi',2,4,-1)
        read(ifiqg ) nqnum, ngpmx ,QpGcut_psi,nqbz,nqi,imx,nqibz
      endif
!! ==== Write, or read past header information, file gwb ====
      ifigwb = fopna('gwb',-1,4)
      rewind ifigwb
      if(job==1) then
        ef0 = 1d99                !dummy
        ifigwn = fopna('gw1',-1,4)
        rewind ifigwn
        ifigwn = fopna('gw2',-1,4)
        rewind ifigwn
        if(rank0) write(ifigwb) nat,nsp,ndima,ndham,alat,qlat,ef0,nqbz,plat,nqnum,nqi
        allocate(ips(nbas),lmxa(nat),bas(3,nat))
        iat = 0
        do  i = 1, nbas
          lmaxa = int( sspec(int(ssite(i)%spec))%lmxa )
          if (lmaxa .gt. -1) then
            iat = iat + 1
            if (iat .gt. nat) call rx('bug in sugw')
            bas(:,iat)=ssite(i)%pos
            lmxa(iat) = lmaxa
          endif
        enddo
        if(rank0) write(ifigwb) lmxa(1:nat), bas(1:3,1:nat)
        deallocate(ips,lmxa)

!!   ... Determine nphimx
        nphimx = 0
        do  i = 1, nbas
          is=ssite(i)%spec
          i_copy_size=size(ssite(i)%pnu)
          call dcopy(i_copy_size,ssite(i)%pnu,1,pnu,1)
          i_copy_size=size(ssite(i)%pz)
          call dcopy(i_copy_size,ssite(i)%pz,1,pnz,1)
          a=sspec(is)%a
          nr=sspec(is)%nr
          z=sspec(is)%z
          rmt(i)=sspec(is)%rmt
          lmaxa = int(sspec(is)%lmxa)
          nmcore=sspec(is)%nmcore
          if (lmaxa .gt. -1) then
            call atwf(0,a,lmaxa,nr,nsp,pnu,pnz,rsml,ehl,rmt(i),z,w,i1,ncore,konfig,ecore,w,w,nmcore)
            nphimx = max(nphimx,i1)
          endif
        enddo
C   ... For band mode
        if (bandmode .and. jobb == 0) then
          allocate(ovvx(ndima,ndima,nsp))
          ovvx = 0d0
        endif

C   ... Atom data (gwa)
        call info(30,1,1,' ... Generate core wave functions (file gwa)', 0,0)
        if(rank0) ifigwa = fopna('gwa',-1,4)
        if(rank0) rewind ifigwa
        do  ib = 1, nbas
          is=ssite(ib)%spec
          i_copy_size=size(ssite(ib)%pnu)
          call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
          i_copy_size=size(ssite(ib)%pz)
          call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
          a=sspec(is)%a
          nr=sspec(is)%nr
          z=sspec(is)%z
          rmt(ib)=sspec(is)%rmt
          lmaxa = int(sspec(is)%lmxa)
          if (lmaxa .gt. -1) then
            call atwf ( 0 , a , lmaxa , nr , nsp , pnu , pnz , rsml , ehl
     .       , rmt ( ib ) , z , ssite(ib)%rv_a_ov0 , i1 , ncore , konfig , ecore , w 
     .       , w ,nmcore)
            allocate(rofi(nr),rwgt(nr),gcore(nr,2,ncore))
            allocate(gval(nr,2,0:lmaxa,nphimx,nsp))
            call dpzero(gval,nr*2*(1+lmaxa)*nphimx*nsp)
!!         Create augmented wave functions for this atom
            call uspecb(0,4,sspec,is,is,lh,rsml,ehl,i)
            call atwf ( 03 , a , lmaxa , nr , nsp , pnu , pnz , rsml , ehl
     .       , rmt ( ib ) , z , ssite(ib)%rv_a_ov0 , nphimx , ncore , konfig , ecore 
     .       , gcore , gval ,nmcore)
!!         Header data for this atom
            b = rmt(ib)/(dexp(a*nr-a)-1d0)
            call radmsh(rmt(ib),a,nr,rofi)
            call radwgt(rmt(ib),a,nr,rwgt)
c --- for band mode
            if (bandmode .and. jobb==0)then
              do  l = 0, lmaxa
                nmx=2
                if (konfig(l) >= 10) nmx=3
                do nn1=1,nmx
                  do nn2=1,nmx
                    do isp=1,nsp
                      call gintsl(gval(1,1,l,nn1,isp),gval(1,1,l,nn2,isp),
     .                a,b,nr,rofi,ovvpp)
                      do mx=1, 2*l+1
                        im1 = nlindx(nn1,l,ipb(ib))+mx
                        im2 = nlindx(nn2,l,ipb(ib))+mx
                        ovvx(im1,im2,isp) = ovvpp
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            endif
            if(rank0)write(ifigwa) z, nr, a, b, rmt(ib), lmaxa, nsp, ncore
            if(rank0)write(ifigwa) konfig(0:lmaxa)
            if(rank0)write(ifigwa) rofi
!!         Write orthonormalized valence wave functions for this atom
            do  l = 0, lmaxa
              do  i = 1, nsp
                if(rank0)write(ifigwa) l,i
                if(rank0)write(ifigwa) gval(1:nr,1,l,1,i)
                if(rank0)write(ifigwa) gval(1:nr,1,l,2,i)
                if (konfig(l) .ge. 10.and.rank0) write(ifigwa) gval(1:nr,1,l,3,i)
              enddo
            enddo
!!         Core wave functions for this atom
            icore = 0
            vshft = vrmt(ib)
!!        As of v6.11, shift is included in v0, passed in vval to
!!         locpot, in routine mkpot.f
            vshft = 0
            do  l = 0, lmaxa
              do  isp = 1, nsp
                do  konf = l+1, mod(konfig(l),10)-1
                  icore = icore+1
                  if(rank0)write(ifigwa) icore, l, isp, konf, ecore(icore)+vshft
                  if(rank0)write(ifigwa) gcore(1:nr,1,icore)
                enddo
              enddo
            enddo
            deallocate(rofi,rwgt,gcore,gval)
          endif
        enddo
        if(rank0)call fclr('gwa',ifigwa)
      else
        call rxi('sugw: bad job',job)
      endif

!! == GW setup loop over k-points ==
      if (lchk .ge. 1 .and. job .eq. 1) then
        ifinormchk = fopna('normchk',-1,0)
        rewind ifinormchk
        write(ifinormchk,849)
  849   format('#     eval          IPW        IPW(diag)    Onsite(tot)   Onsite(phi)      Total')
      endif

!! --- Band mode --- not tested recently.
      if (bandmode .and. jobb .eq. 0) then
        if (pwmode .ne. 0) call rx('band mode not ready for APW')
        ndimh = ldim + napw
        allocate(iiyf(ndimh),ibidx(ndimh,nsp),ovv(ndimh,ndimh),
     .    cphi_p(ndima,ndimh,nsp),evl_p(ndham,nsp))
C   ... Open LBAND files
        do  isp = 1, nsp
          if (isp.eq.1) iflband(isp) = fopnx('LBAND.UP',2,6,-1)
          if (isp.eq.2) iflband(isp) = fopnx('LBAND.DN',2,6,-1)
          write(iflband(isp)) ndimh,nqnum
          write(iflband(isp)) plat,qlat
        enddo
C   ... Read SYML file
        nlinemax = 50
        allocate(nqq(nlinemax),qq1(1:3,nlinemax),qq2(1:3,nlinemax))
        ifsyml = fopnx('SYML',2,1,-1)
        nline = 0
        do
          nline = nline + 1
          read(ifsyml,*,err=601,end=601) nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
        enddo
  601   continue
        call fclose(ifsyml)
        nline = nline - 1
        qpos  = 0d0
        iline = 1
        if (nline .eq. 0) call rx('sugw: no lines in SYML file')
      elseif (bandmode .and. jobb==5) then
        ifqeigen = fopnx('Qeigval',2,1,-1)
        read(ifqeigen,*) nqnum
        ifeigen = fopnx('eigval',2,4,-1)
        if(rank0) write(ifeigen) ndham,nqnum,nsp
      endif

!! Write, or read ngpmx
      if(rank0) then
         ifigwb = fopna('gwb',-1,4)
         if (job .eq. 1) then
         write(ifigwb) ngpmx
         else
         read(ifigwb)
         endif
      endif

!! mar2012 rotwv mode. on regular mesh, only eigenfunctions for nqibz are required.
      nqtot=nqnum
      nnn=nqibz

!! --- Evecs and matrix elements of vxc for irr qp ---
!!    Note: this routine should use only irr qp.
      if (lsig > 0) then
        if (job==1) then
        if (.not. cmdopt('--novxc',7,0,strn)) then
           ifievec = fopna('evec',-1,4)
           rewind ifievec
           if(rank0) write(ifievec) ndham, nsp, nnn !, nqnum
           ifiv = fopna('vxc',-1,4)
           rewind ifiv
           if(rank0) write(ifiv) ndham, nsp, nnn
        endif
        endif
      endif

!! == Main loop for eigenfunction generation ==
      if (ipr .ge. 20) then
         if (job .eq. 1 .and. bandmode) call info2(20,1,1,
     .        ' ... Make LDA eigenvalues, %i qp (file eigval)', nqtot,0)
         if (job .eq. 1 .and. .not. bandmode) call info2(20,1,1,
     .        ' ... Make LDA w.f. and matrix elements, %i qp: file(s) gwb'//
     .        '%?#n#, vxc, evec##',nqtot,lsig)
         if(sham%scaledsigma/=1d0.and.sigmamode) then
           write(6,*)' Scaled Sigma method: ScaledSigma=',sham%scaledsigma
         endif   
      endif
      iqibz=0

!!  === job divider apr2013takao ===
      nbalance = nqi*nsp/numproc 
      if(mod(nqi*nsp,numproc)/=0) nbalance = nbalance+1
      if(iprint()>10) print *,'nbalance=',nbalance
      iproc = numproc-1         !from bottom
      allocate(iprocq(nqi,nsp))
      iqq  = 0
      iadd = 0
      do iq = 1, nqtot
         read(ifiqg)  q,  ngp, irr   
         read(ifiqg)  
         if(irr==0) cycle
         iqq = iqq + 1
         do isp  = 1, nsp
            iadd = iadd+1
            iprocq(iqq,isp)=iproc
            if(iadd == nbalance) then
               iproc = iproc-1
               iadd = 0
            endif  
         enddo
      enddo
      nrr=iqq
      if(nqi/=nrr) call rxi('sugw: nqi/=nrr=',nrr)

      if(rank0) then
         ifiproc = fopnx('lmfgw_kdivider',2,2,-1)
         call fextg(ext)
         write(ifiproc,"(a,'    ! extension ')") trim(ext)
         write(ifiproc,"(3i10,'    ! nqi nsp numproc')") nqi, nsp, numproc
         do isp=1,nsp
         do iqq=1,nqi
            write(ifiproc,"(3i9,'  ! iqq isp irank')") iqq, isp, iprocq(iqq,isp) 
         enddo
         enddo
         call fclose(ifiproc)
      endif    
!!
      rewind ifiqg
      read(ifiqg) nqnum, ngpmx ,QpGcut_psi,nqbz,nqi,imx,nqibz
      if(pwmode<10) call shortn3_initialize(qlat,rlatp,xmx2)
      if(debug)print *,' goto 1001 procid=',procid
      iqq=0
      do 1001 iq = 1, nqtot
        read (ifiqg)  q,  ngp, irr   ! q, and number of G vectors for
        ! eigenfunction expansion at this q
        ! expansion of the Coulomb interaction
        if(irr==1) iqibz=iqibz+1
        lwvxc = lsig>0 .and. job==1 .and. (.not.bandmode).and. iqibz<=nqibz
        if (cmdopt('--novxc',7,0,strn)) lwvxc = .false.
        if (jobb .eq. 5) then
          read(ifqeigen,*) q
          goto 1021
        endif
        if(allocated(ngvecp)) deallocate(ngvecp)
        allocate(ngvecp(3,max(ngp,1))) 
        read (ifiqg)  ngvecp
        if (job .ne. 1) deallocate(ngvecp) 
!! calculate only for irr=1 case
        if(irr==0) then
          cycle
        endif
 1021   continue

!! cycle when this iqq is not for procid apr2013takao
        iqq = iqq+1
        if(iprocq(iqq,1)/=procid .and. iprocq(iqq,nsp)/=procid) then
          if(allocated(ngvecp)) deallocate(ngvecp)
          cycle
        endif  
!!
        qp = q
cnnn
        if(pwmode>0.and.pwmode<10) then
        ppin=matmul(transpose(plat),qp) !basis on the qlat coordinate. qp in Cartesian.
        call shortn3(rlatp,xmx2,ppin,noutmx, nout,nlatout)
        do iout=1,nout
           write(*,"(a,3i5,f10.4,3f8.4)")'sugw:rrrrn1 =',nlatout(:,iout),
     &         sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2),
     &         matmul(qlat(:,:),ppin+nlatout(:,iout))
        enddo
        endif

c        write(stdo,"(' procid:  iq=',2i4,'  q=',3f7.3)") procid, iq,qp
!! keep comment for noshorbz case.
c$$$        if(.not.noshorbz()) then
c$$$          call shorbz(q,qp,qlat,plat)
c$$$        else
c$$$           qp=q
c$$$        endif
c$$$        write(stdo,"(' iq=',i4,'  q=',3f7.3,'  shortened to',3f7.3)") iq,q,qp

!!!  ... For this qp, G vectors for PW basis and hamiltonian dimension
        if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then
          pwgmin = dsqrt(pwemin)
          pwgmax = dsqrt(pwemax)
          call pshpr(1)
          qqq=0d0 !call dpzero(xx,3)
          if (mod(pwmode/10,10) == 1) qqq=qp !call dpcopy(qp,xx,1,3,1d0)
          call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,0,napw,dum,dum,dum,dum)
          call poppr
          allocate(igv2(3,napw))
          allocate(kv_iv(3*napw))
          call pshpr(iprint()-10)
c          call gvlst2(alat, plat, xx, 0, 0, 0, pwgmin, pwgmax, 0, 2, 0, napw, kv_iv, dum, dum, igv2)
          call gvlst2(alat, plat, qqq, 0, 0, 0, pwgmin, pwgmax, 0, 2, napw, napw, kv_iv, dum, dum, igv2)
cnnn
          if (pwmode<10) then
             do iapw=1,napw
                igv2(:,iapw)=igv2(:,iapw)+nlatout(:,1)
             enddo
          endif

          if (allocated(kv_iv)) deallocate(kv_iv)
          call poppr
          ndimh = ldim + napw
          if (mod(pwmode,10) .eq. 2) ndimh = napw
          if (ndimh .gt. ndham) then
            call fexit2(-1,111,'%N Exit -1 : BNDFP: '//
     .      'ndimh=%i exceeds ndham=%i.  Set '//
     .      'HAM_NPWPAD=(ndimh-ndham) or larger.',ndimh,ndham)
          endif
        else
          ndimh = ldim
          allocate(igv2(1,1)) !dummy
        endif

!! this is old comment --->   ... Must pass shortened qp to pwmat in APW case
c$$$        if (napw .gt. 0 .and. loldpw .lt. 2) then
c$$$          loldpw = loldpw + 2
c$$$        endif
        allocate(ham(ndimh,ndimh),ovl(ndimh,ndimh),evec(ndimh,ndimh))
        allocate(vxc(ndimh,ndimh))
        allocate(cphi(ndima,ndimh,nsp),cphin(2,ndimh,nsp))
        if (bandmode .and. jobb .eq. 0) then
          allocate(geig(max(ngp,1),ndimh,nsp))
        endif
        if(debug) print *,' goto hambl loop'
        do 1002 isp = 1, nsp
!! ==== cycle when isp for iqq is not for procid apr2013takao ===
          if(iprocq(iqq,isp)/=procid) cycle
!!   --- vxc <- LDA Hamiltonian without vxc for this qp ---
          alfa = 0
          call hambl ( 0 , nbas , ssite , sspec , slat , sham , isp , qp 
     .       , k1 , k2 , k3 , spotx , vconst , sv_p_osigx , sv_p_otaux , sv_p_oppix 
     .       , lcplxp , alfa , ndimh , napw , igv2 , vxc , ovl , w )
          if(debug) print *,' goto hambl vxc '
!!  --- LDA Hamiltonian and overlap matrices for this qp ---
          call hambl ( 0 , nbas , ssite , sspec , slat , sham , isp , qp 
     .       , k1 , k2 , k3 , smpot , vconst , sv_p_osig , sv_p_otau , sv_p_oppi 
     .       , lcplxp , alfa , ndimh , napw , igv2 , ham , ovl , w )
          if(debug) print *,' goto hambl lda o '
!!   --- vxc <- LDA vxc ---
          vxc = ham - vxc
!!   --- Write matrix elments of vxc to disk (Faleev May 2002) ---
          if (lwvxc) then
c            write(6,"('lwvxc=T write vxc: iq q=',i3,3f8.4)")iq, q
            write(ifiv) ndimh, ldim,qp !Add ldim sep2012
            write(ifiv) q,vxc
          endif
!!   --- LDA + sigma Hamiltonian for this qp ---



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          if(sigmamode) then
             allocate(sene(ndimsig,ndimsig))
             i = lrsig*10 + 1
             hreal = 0
             if (mod(i/10,10) >= 4) hreal = 1 !this is lrsig in hambles.
             mode = 100000 + 4000 + 40*(1-hreal) + 000 + 10 !hreal=1 -> sig(T) is real
             if(debug) print *,'AAAAA2 goto bloch: i q nl plat mxorb=',lrsig,mode,qp,ndimsig,ival ( sham%iv_a_ontabs , nbas + 1 )

             call bloch (mode,qp,globalvariables%nl,plat,globalvariables%mxorb,sham%iv_a_oindxo,1
     &            , ival ( sham%iv_a_ontabs , nbas + 1 )
     &            , sham%iv_a_oiaxs , sham%rv_a_ohrs , sham%ndhrs , isp , nsp , ndimsig , ndimsig 
     &            , 0 , ndimsig, 0 , ndimsig , 0 , sene, iwdummy, iwdummy )
!!     case2 not work well yet. We need think about pair tables as in bloch.
c            allocate(sene(ndimsig,ndimsig))
c     call seneinterp2(qp,nbas,nk1,nk2,nk3,plat, !npair,nlat,qwgt,npairmx,qsmesh,
c     i        sfz(:,:,:,:,:,isp),ldim,ndimsig, 
c     o        sene) 
             if(debug) print *,'end of bloch sene'
          endif

          lpdiag = 0
          jx = 0
          call hambl ( jx , nbas , ssite , sspec , slat , sham , isp ,
     .         qp , k1 , k2 , k3 , smpot , vconst , sv_p_osig , sv_p_otau 
     .         , sv_p_oppi , lcplxp , alfa , ndimh , napw , igv2, ham, ovl, ham)

          if(sigmamode) then
             allocate( ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
             ovlmtoi = ovl(1:ndimsig,1:ndimsig)
             call matcinv(ndimsig,ovlmtoi)
             ovliovl = matmul(ovlmtoi,ovl(1:ndimsig,1:ndimh))
             deallocate(ovlmtoi)
             allocate(senex(ndimh,ndimh))
             senex = matmul(transpose(dconjg(ovliovl)), matmul(sene,ovliovl))
             deallocate(ovliovl,sene)
             ham(:,:) = ham(:,:) + sham%scaledsigma * senex
             deallocate(senex)
          endif
c$$$          if(sigmamode.and.newsigmasw()) then
c$$$            allocate(sene(ndimh,ndimh))
c$$$            call seneinterp2(qp,nbas,nk1,nk2,nk3,plat,
c$$$     i       sfz(:,:,:,:,:,isp),ldim,ndimh, 
c$$$     o       sene)
c$$$            ham=ham+sene
c$$$            deallocate(sene)
c$$$          else
c$$$            ispc = min(isp,nspc)
c$$$            isize_ham=(/size(ham,dim=1),size(ham,dim=2),1/)
c$$$            call hambls ( i , nbas , ssite , sspec , slat , sham , isp ,
c$$$     .        ispc , qp , k1 , k2 , k3 , sham%rv_p_oqsig , nqsig , smpot , vconst 
c$$$     .        , sv_p_osig , sv_p_otau , sv_p_oppi , lcplxp , 0 , alfa , ndimh 
c$$$     .        , napw , igv2 , ham , ovl , isize_ham, i2, ndimsig )
c$$$          endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



          if (i .eq. -1) lpdiag = 2
          if(debug) print *,' end of generation of Hamiltonian'

!!   --- Branch job = 1 : make cphi, matrix elements ---
          if(job .eq. 1) then !   --- Diagonalize ---
            if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
            write(6,"(' sugw:  kpt isp=',i8,i2,' of ',i8, ' k= ',3f9.5, ' ndimh= ',i5,
     &      ' irank=',i4, ' lwvxc=',l,' lpdiag=',i1)")  iq,isp,nqtot,qp,ndimh,procid,lwvxc,lpdiag
c            if(oncewrite(5)) print *,'vvvvv lpdiag=',lpdiag
            if (lpdiag .eq. 0) then
              if (nspc .eq. 2) then
                call rx('diagonalization not ready for nspc=2')
              endif
              allocate(ww_rv(11*ndimh))
              epsovl = sham%oveps
              if (epsovl .le. 0) then
                call zhev ( ndimh , ham , ovl , .true. , .true. , ndimh , 1d60 
     .          , nev , ww_rv , .false. , - 1 , evl ( 1 , isp ) , evec )
              else
                nevl = -1
                call dvset ( ww_rv , 1 , 1 , 99999d0 )
                call zhevo ( ndimh , ndimh , ham , ovl , .true. , ndimh , 1d60 
     .          , epsovl , nevl , nev , evl ( 1 , isp ) , ww_rv , ndimh , 
     .          evec )
                eomin = dval ( ww_rv , 1 )
                call info5(30,0,0, ' Overlap''s smallest eigenvalue: %;3g.  '//
     .          '%?#(n>0)#H dim reduced from %i to %i#H dim not reduced#',
     .          eomin,ndimh-nevl,ndimh,nevl,0)
              endif
              deallocate(ww_rv)
              if(debug) print *,'wwwww ndimh nev=',ndimh,nev
              call prtev(evec,ndimh,evl(1,isp),ndimh,1d60,nev)
            elseif (lpdiag .eq. 2) then
              evec = ovl
#ifdef USEOPTIONAL
              call phmbls(mode=2,ndimh=ndimh,eval=evl(1,isp),zhz=ham)
#else
              call phmbls(2,ndimh,evl(1,isp),w,w,w,w,w,ham)
#endif
              nev = ndimh
              call prtev(evec,ndimh,evl(1,isp),ndimh,9d9,nev)
            else
              call rxi('sugw not ready for lpdiag=',lpdiag)
            endif
C       Pad evals between ndimh and ndham with a large positive number
C       to avoid mixing up integration routines
            if (ndham .gt. nev) call dvset(evl(1,isp),1+nev,ndham,99999d0)
            if (bandmode .and. jobb .eq. 5) then
              write(ifeigen) evl(1:ndham,isp)
              if (mod(iq,10) .ne. 1) call poppr
              cycle
            endif
            if (mod(iq,10) .ne. 1) call poppr
            if(debug) write(6,"(' sugw:procid iq isp lwvxc= ',3i3,' ',l)")procid, iq,isp,lwvxc 
            if (lwvxc) then
              ifievec = fopna('evec',-1,4)
              write(ifievec) q, evec(1:ndimh,1:ndimh)
            endif

            if(magexist) then
               if(isp==1) evl(1:ndimh,isp)=evl(1:ndimh,isp) - vnow/2d0
               if(isp==2) evl(1:ndimh,isp)=evl(1:ndimh,isp) + vnow/2d0
            endif   

            if(debug) write(6,"(' ppppppp 111111111')")
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$!! test for takao wave function rotater si 5x5x5.
c$$$            if(.false.) then
c$$$c       if(iq==8) then
c$$$              ifievv=1012
c$$$              open( ifievv, file='evecx' )
c$$$              ifiqbz=1025
c$$$              open(ifiqbz,file='QBZ')
c$$$              read(ifiqbz,*) nqzz
c$$$              allocate(qzz(3,nqzz))
c$$$              do i=1,nqzz
c$$$                read(ifiqbz,*) qzz(:,i)
c$$$              enddo
c$$$
c$$$              allocate(evecout(ndimh,ndimh))
c$$$              do iqzz = 1,nqzz
c$$$                qtarget = qzz(:,iqzz)
c$$$c
c$$$c                call rotwv(q,qtarget,ndimh,napw,ndimh, plat,qlat,evec,evecout,ierr)
c$$$                call rotwv(q,qtarget,ndimh,napw,ndimh, evec,evecout,ierr)
c$$$                if(ierr/=0) cycle
c$$$
c$$$                do
c$$$                  read(ifievv,*,end=1019) iiiii
c$$$                  read(ifievv,*) qrr, ndimh
c$$$                  print *,'qrr ndimh=',qrr,ndimh
c$$$                  allocate(evecr(ndimh,ndimh))
c$$$                  print *,'qrr ndimh xxx=',qrr,ndimh
c$$$                  do j= 1,ndimh
c$$$                    do i= 1,ndimh
c$$$                      read(ifievv,*) axx,bxx
c$$$                      evecr(i,j)=dcmplx(axx,bxx)
c$$$                    enddo
c$$$                  enddo
c$$$                  if(sum(abs(qrr-qtarget))<1d-8) exit
c$$$                  deallocate(evecr)
c$$$                enddo
c$$$                rewind ifievv
c$$$c$$$                if(.not.noshorbz()) then
c$$$c$$$                  call shorbz(qtarget,qxxx,qlat,plat)
c$$$c$$$                else
c$$$c$$$                  qxxx=qtarget
c$$$c$$$                endif
c$$$                qxxx=qtarget
c$$$                write(1013,"(i10)") 11111
c$$$                write(1013,"(3f8.3,i10)") qtarget,ndimh
c$$$c           write(1013,"(i10,3f8.3,i10)") 11111,qout,ndimh
c$$$                do j=1,ndimh
c$$$                  do i=1,ndimh
c$$$                    if(abs(evecout(i,j))+abs(evecr(i,j))>1d-4) then
c$$$                      write(1013,"(2i4,2d13.5,2x,2d13.5,2x,2d12.4,2x,2d12.4,3f8.3)") 
c$$$     &                i,j,evecout(i,j), evecr(i,j), evecout(i,j)/evecr(i,j),
c$$$     &                abs(evecout(i,j)), abs(evecr(i,j)),qxxx
c$$$                    endif
c$$$                  enddo
c$$$                  write(1013,*)
c$$$                enddo
c$$$
c$$$                deallocate(evecr)
c$$$              enddo
c$$$              deallocate(evecout)
c$$$              stop 'test end xxxxxxxx'
c$$$ 1019         continue
c$$$              stop 'uuuuuuuuuuuuuuu'
c$$$            endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$!! takao write eigenfunctions to file fort.ifievv.
c$$$            if(.false.) then
c$$$              write(ifievv,"(i10)") 11111
c$$$              write(ifievv,"(3f8.3,i10)") q,ndimh
c$$$              do j=1,ndimh
c$$$                do i=1,ndimh
c$$$                  write(ifievv,"(2i4,2d13.5,2x,d13.5)") i,j,evec(i,j),abs(evec(i,j))
c$$$c          write(ifievv,"(2d23.15)") evec(i,j)
c$$$                enddo
c$$$              enddo
c$$$            endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!!   --- Project wf into augmentation spheres, Kotani conventions ---
            nlmax = globalvariables%mxorb / globalvariables%nkaph
            allocate(aus_zv(abs(-nlmax*ndham*3*nsp*nbas)))
            if (-nlmax*ndham*3*nsp*nbas<0) aus_zv(:)=0.0d0
            if(debug) write(6,"(' ppppppp goto makusq')")
            call makusq ( 1 , ssite , sspec , slat , sham , nbas , nbas , 
     .      0 , nlmax , ndham , ndimh , napw , igv2 , nev , nsp , 
     .      nspc , isp , 1 , qp , evec , ppn , aus_zv )
            call gwcphi ( ssite , sspec , isp , nsp , nlmax , ndham , nev 
     .      , nbas , ipb , lmxax , nlindx , ndima , ppn , aus_zv , cphi 
     .      ( 1 , 1 , isp ) , cphin ( 1 , 1 , isp ) )
            if (allocated(aus_zv)) deallocate(aus_zv)
            if(debug) write(6,"(' ppppppp end of gwcphi')")

!! We keep note in the followings, but be careful (may contain bugs)...
!!  --- Overlap of IPWs, PW expansion of eigenfunctions pwz ---
!!      The IPW basis consisting of PWs with holes knocked out of spheres,
!!      IPWs must be orthogonalized: IPW -> IPWbar; see PRB 76, 165106.
!!         |IPWbar_G> = sum_G2 |IPW_G2> O^-1_G2,G1, where O_G1,G2=<IPW_G1|IPW_G2>
!! == Definitions ==
C     * ppovl = overlap of IPWs (Generated by pwmat and pwmat2)
C         ppovl_G1,G2 = O_G1,G2 = <IPW_G1 IPW_G2>
C     * pwh = PW expansion of basis function (Generated by pwmat2 only)
C         basis_j> = sum_G2 PWH_G2,j |IPW_G2>
C     * Matrix elements (overlap) of IPW and basis function (Generated by pwmat only)
C         phovl_G1,j = sum_G2 ppovl_G1,G2 pwh'_G2,j  (matrix form PHOVL = O * PWH')
C         Note: phovl is only used as an intermediate construction, old branch
C     * Note: pwh' is expanded to the LMTO cutoff gmax while
C             pwh  is expanded to the GW cutoff QpGcut_psi
C       Thus O^-1 PHOVL will not identically recover PWH.
C       The original branch (loldpw=0) uses effectively PWH'; the new one uses PWH.
C       This is a major distinction between the two (see Remarks)
C     * PW expansion of eigenfunction:
C         |psi_n> = sum_j z_jn |basis_j>
C     * Define pwz_G,n = sum_j PWH_G2,j z_jn  (in matrix form: PWZ = PWH Z)
C       Then
C         |psi_n> = sum_j,G2 z_jn PWH_G2,j |IPW_G2> = sum_G2 PWZ_G2,n |IPW_G2>
C       Overlap of IPW and eigenfunction:
C         PZOVL_G1,n = <IPW(G1) psi_n> = sum_G2 O_G1_G2 PWZ_G2,n
C         PZOVL = O * PWZ (matrix form) <--- old
C         Note: pzovl is only used as an intermediate construction, old branch
            if(debug) print *,'ppppppppppp ngp=',ngp
            if (ngp .gt. 0) then
              allocate(ppovl(ngp,ngp),pwz(ngp,ndimh))
              allocate(phovl(ngp,ndimh))
C         Pass qx to pwmat (or pwmat2):
C         qx = (unshortened) q if 2s digit loldpw = 0
C         qx = (shortened)  qp if 2s digit loldpw = 1
C         Old convention: call pwmat
c          if (mod(loldpw,2) .eq. 0) then

!!  We have  q+G(igvx; internal in pwmat) = qp + G(igv2)  
!!  Thus, we have
!!           igv(internally in pwmat) = igv2 + qlatinv*(qp-q)
!!            inn = qlatinv*(qp-q)
!!!! NOTE: from 20Sep2012, we set qp=q (shorbz is in hambl.F. In other words, q supplied to hambl
!!!! do not need to be short enough).
              if(sum(abs(qp-q))>1d-8) stop 'sugw:qp/=q; qp=p sep2012'
              inn=0
              if(debug) print *,'goto pwmat'
              call pwmat ( slat , ssite , sspec , nbas , ndimh , napw,
     .          igv2, sham%iv_a_oindxo , q , ngp , nlmax , ngvecp , gmax , inn, 
     .          ppovl, phovl )
              if(debug) print *,'sss: ppovl=',sum(abs(ppovl))
              if(debug) print *,'sss: phovl=',sum(abs(phovl))
              if(debug) print *,'sss: evec =',sum(abs(evec))
              if(debug) print *,'sss:       '
              call zgemm('N','N',ngp,ndimh,ndimh,(1d0,0d0),phovl,ngp, evec,ndimh,(0d0,0d0),pwz,ngp)
              if(debug) print *,'sss: pwz =',sum(abs(pwz))
              if(debug) print *,'sss:       '
              deallocate(phovl)
              if (lchk .ge. 1) then
                allocate(pzovl(ngp,ndimh))
                pzovl = pwz
                allocate(ppovld(ngp)) ! extract diagonal before ppovl overwritten
                do  i = 1, ngp
                  ppovld(i) = ppovl(i,i)
                enddo
              endif
              call zqinvb('h',ppovl,ngp,ngp,ndimh, pwz,ngp, ierr)
              if (ierr .ne. 0) call rx('zqinvb failed to invert ppovl')
            endif
            if(debug) print *,' goto gwb write...'
C   --- File output, job=1 ---
            ifigwb = fopna('gwb',-1,4)
            write (ifigwb) q, ndimh, iq
            if(debug) write (6,"('q ndimh=',3f10.5,i10)") q, ndimh
            write (ifigwb) evl(1:ndimh,isp),cphi(:,:,isp)
            write (ifigwb) ngp !takao feb2012 ,ngc
            write (ifigwb) ngvecp, pwz !takao feb2012
            ifigw1 = fopna('gw1',-1,4)
            write(ifigw1) q, (evl(i,isp),i=1,ndimh)
C       call prmx('e(H+sigma-vxc)',evl,ndimh,ndimh,1)

C   ... Overlap checking.   Define:
C       Interstitial part of eigenfunction overlap:
C         <psi_n||psi_n'>
C         = sum_G1,G2 (pwz_G1,n|IPW_G1>)+  (pwz_G2,n'|IPW_G2>)
C         = sum_G1,G2 (pwz_G1,n)+ ppovl_G1,G2 (PWZ_G2,n')
C         = (PWZ)+ O (PWZ) = (PZOVL)+ (PWZ)  (old style)
            if (lchk .ge. 1 .and. ngp .gt. 0) then
              allocate(testc(ndimh,ndimh),testcd(ndimh))
c$$$          if (mod(loldpw,2) .eq. 1) then
c$$$            allocate(pzovl(ngp,ndimh))
c$$$            call zgemm('N','N',ngp,ndimh,ngp,dcmplx(1d0,0d0),
c$$$     .        ppovl,ngp,pwz,ngp,dcmplx(0d0,0d0),pzovl,ngp)
c$$$          endif
              call zgemm('C','N',ndimh,ndimh,ngp,(1d0,0d0),
     .        pzovl,ngp,pwz,ngp,(0d0,0d0),testc,ndimh)
              deallocate(pzovl)
C          call zprm('(PWZ)+^-1 O_i^-1 (PWZ)',2,testc,ndimh,ndimh, ndimh)
              do  i = 1, ndimh
                testcd(i) = sum(dconjg(pwz(:,i))*ppovld*pwz(:,i))
              enddo
              deallocate(ppovld)
C         xx(1) = sum over all augmentation w.f.  cphi+ ovl cphi
C         xx(2) = sum over augmentation phi only.
C         xx(3) = IPW contribution to phi+ phi
C         xx(4) = IPW contribution to phi+ phi, using diagonal part only
              ifinormchk = fopna('normchk',-1,0)
              if (abs(sum(q-qp)) .gt. 1d-10) then
                write(ifinormchk,555) iq,q,qp
              else
                write(ifinormchk,555) iq,q
              endif
  555         format('# iq',i5,'   q',3f12.6:'  shortened q',3f12.6)
              do  i1 = 1, ndimh
                xx(1) = cphin(1,i1,isp)
                xx(2) = cphin(2,i1,isp)
                do  i2 = 1, ndimh
C         xx(1) = sum(dconjg(cphi(1:ndima,i1,isp))*cphi(1:ndima,i2,isp))
C         xx(2) = sum(dconjg(cphi(1:nchan,i1,isp))*cphi(1:nchan,i2,isp))
                  xx(3) = testc(i1,i2)
                  xx(4) = testcd(i1)
                  if (i1 .eq. i2) then
                    write(ifinormchk,'(f12.5,5f14.6)')
     .              evl(i1,isp),xx(3),xx(4),xx(1),xx(2),xx(1)+xx(3)
                  endif
                enddo
              enddo
              deallocate(testc,testcd)
              write(ifinormchk,*)

C   ... Band mode
              if (bandmode) then
                geig(:,:,isp)= pwz
                if (isp .eq. 1) then
                  if (iq .eq. 1) q_p=q
                  qpos = qpos + sqrt( sum((q-q_p)**2) )
                  if (endofline) then
                    iline = iline + 1
                    endofline=.false.
                  endif
                  if (sum(abs(q-qq2(:,iline)))<1d-6) endofline = .true.
                endif

                if (iq .eq. 1) then
                  do  j = 1, ndimh
                    ibidx(j,isp) = j
                  enddo
C           Should not be necessary, but avoids a bug in the DEC alpha compiler
                  ngp_p = ngp
                else
                  call tcn('bndconn')
                  call bndconn_v2( alat,plat,qlat,
     i            nbas, rmt, bas, ndimh, ndima,
     i            evl(1,isp),  ngp, ngp_p, ngvecp,ngvecp_p,
     i            geig(1,1,isp),   geig_p(1,1,isp),
     i            cphi(1,1,isp), cphi_p(1,1,isp), ovvx(1,1,isp),
     o            iiyf,ovv)
C           Continuous band index:
                  ibidx(1:ndimh,isp) = iiyf(ibidx(1:ndimh,isp))
                  write(96+isp,'(" ",    255i4)') iq,(j,      j=1,ndimh)
                  write(96+isp,'("     ",255i4)') (iiyf(j),j=1,ndimh)
                  write(96+isp,'("     ",255i4)')
     .            (int(10*ovv(iiyf(iix),iix)),iix=1,ndimh)
                  call tcx('bndconn')
                endif
                write(995,"(i2,3f9.5,' ',f9.5,' ',255(i3,f10.5))") iline,
     .          q(1:3), qpos,
     .          (ibidx(ibr,isp),evl(ibidx(ibr,isp),isp), ibr=1,ndimh)
                write(iflband(isp)) iline, q(1:3), qpos
     .          ,(ibidx(ibr,isp),evl(ibidx(ibr,isp),isp),ibr=1,ndimh)

!!     --- Save eigenfunctions of previous loop ---
                if(iq/=1 .and. isp.eq.nsp) deallocate(geig_p, ngvecp_p)
                if (isp .eq. nsp) then
                  allocate(geig_p(ngp,ndimh,nsp), ngvecp_p(3,ngp))
                  cphi_p = cphi
                  geig_p = geig
                  evl_p  = evl
                  ngvecp_p = ngvecp
                  ngp_p    = ngp
                  q_p      = q
                endif
              endif                   !------------- end of band mode
            endif                   !------------- end of lchk=1
            if (ngp .gt. 0) deallocate(ppovl,pwz)
            if (isp .eq. nsp) deallocate(ngvecp)!takao feb2012 ,ngvecc)
!!   ... Subtract <psi|vxc|psi> from starting evl=e(lda+sigma-vxc)
!!       evl should contain e(lda+sigma-vxc)
            allocate(testc(ndimh,ndimh))
            call zgemm('C','N',ndimh,ndimh,ndimh,(1d0,0d0),
     .      evec,ndimh,vxc,ndimh,(0d0,0d0),testc,ndimh)
            do i1 = 1, ndimh
            do i2 = 1, ndimh
              evl(i1,isp) = evl(i1,isp) - testc(i1,i2) * evec(i2,i1)
            end do
            end do
C       call prmx('e(H+sigma-vxc-vxc)',evl,ndimh,ndimh,1)
C   ... Save to file gw2 energies e(lda+sigma-vxc-vxc)
            ifigw2 = fopna('gw2',-1,4)
            write(ifigw2) q, (evl(i,isp),i=1,ndimh)
            deallocate(testc)
c$$$C   --- Branch job = 2 : matrix elements of ham without vxc ---
c$$$          else
c$$$
c$$$C   ... Diagonal matrix elements of <psi|H(lda)-vxc|psi>
c$$$            allocate(testc(ndimh,ndimh))
c$$$C       Read eigenvectors from gwb
c$$$            ifi = fopna('gwb',-1,4)
c$$$            read (ifi)
c$$$            read (ifi) evec
c$$$            read (ifi)
c$$$            read (ifi)
c$$$C       evl = <psi|H(no XC)|psi>
c$$$            call zgemm('C','N',ndimh,ndimh,ndimh,(1d0,0d0),
c$$$     .      evec,ndimh,ham,ndimh,(0d0,0d0),testc,ndimh)
c$$$            call dpzero(evl,ndimh)
c$$$            do  i1 = 1, ndimh
c$$$              do  i2 = 1, ndimh
c$$$                evl(i1,1) = evl(i1,1) + testc(i1,i2) * evec(i2,i1)
c$$$              end do
c$$$            end do
c$$$C       call prmx('e(H-vxc)',evl,ndimh,ndimh,1)
c$$$C       File write
c$$$            ifi = fopna('gw2',-1,4)
c$$$            write(ifi) q, (evl(i,1),i=1,ndimh)
c$$$            deallocate(testc)
c$$$
          endif
 1002  continue             ! Loop over spins
        deallocate(ham,ovl,evec,vxc,cphi,cphin,igv2)
        if (bandmode .and. jobb .eq. 0) deallocate(geig)
 1001 continue ! ===== Loop over qp ===============================

      if (jobb .ne. 5) call fclr(' ',ifiqg)
      if (.not. bandmode) then
        ifigwb = fopna('gwb',-1,4)
        call fclose(ifigwb)
      endif
      deallocate(evl,nlindx)
      if (bandmode) then
        if (jobb .eq. 5) call rx0('lmfgw mode=5 (Qeigval--->eigval)')
        if (job  .eq. 1) call rx0('lmfgw mode=4 (bandmode)')
      endif
      call fclr('gw1',fopna('gw1',-1,4))
      call fclr('gw2',fopna('gw2',-1,4))
      call fclr('gwb',fopna('gwb',-1,4))
      if (.not. cmdopt('--novxc',7,0,strn)) then
        if (lsig .gt. 0) then
          call fclr('vxc',fopna('vxc',-1,4))
          call fclr('evec',fopna('evec',-1,4))
        endif
      endif
      call tcx('sugw')
      end subroutine sugw



      subroutine ioaindx(npqn,lmxax,nbas,ndima,nlindx,ifi)
C-  File I/O of NlAindx
C ----------------------------------------------------------------------
Ci Inputs
Ci   npqn  :leading dimension of nlindx
Ci   lmxax :second dimension of nlindx
Ci   nbas  :size of basis
Ci   ndima :number of augmentation channels
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Co Inputs/Outputs
Cio  nlindx:pointer to augmentation channels by site
Cio        :(ifi>0) input
Cio        :(ifi<0) output
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npqn,lmxax,nbas,ndima,ifi,nlindx(npqn,0:lmxax,nbas)
C ... Local parameters
      character outs*80
      integer i,ipqn,l,ib,ii

C --- File read ---
      if (ifi .gt. 0) then
        nlindx = -1
        rewind ifi
        read(ifi,'(a)') outs
        read(ifi,*) i
C       If passed ndima>0, check that it matches file value
        if (ndima .gt. 0 .and. i .ne. ndima)
     .  call rx('ioaindx: file mismatch NLAindx')
        ndima = i
        do  i = 1, ndima
          read(ifi,'(a)',err=101,end=101) outs
          read(outs,*) ipqn,l,ib,ii
          nlindx(ipqn,l,ib) = ii
        enddo
  101   continue

C --- File write ---
      else
        call rx('ioaindx: file write not implemented')
      endif

      end subroutine ioaindx




      subroutine wlattc(alat,plat,nbas,nat,ipb,lmxax,lmxa,nsp,konf)

C- Write LATTC file
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nbas  :size of basis, including floating orbitals sites
Ci   nat   :size of true basis (exclude floating orbitals sites)
Ci   ipb   :index to true basis (excluding floating orbitals)
Ci         :given site index including those orbitals
Ci   lmxax :global maximum of lmxa
Ci   lmxa  :augmentation l-cutoff
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   konf  :principal quantum numbers
Co Outputs
Co   file LATTC is written to disk
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   05 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nat,nsp,lmxax,lmxa(nbas),konf(0:lmxax,nbas,nsp)
      integer ipb(nbas)
      real(8) :: plat(3,3),alat
C ... Local parameters
      integer ifi,ib,isp,fopnx
C     print *, ' --- writing LATTC ---'
C     don't use this ... non portable, and can fail on some machines
C     ifi=6661
C     open (ifi, file='LATTC')
      ifi = fopnx('LATTC',2,2,-1)

      write(ifi,"(e24.16)") alat
      write(ifi,"(3e24.16)") plat(1:3,1)
      write(ifi,"(3e24.16)") plat(1:3,2)
      write(ifi,"(3e24.16)") plat(1:3,3)
      write(ifi,*) ' -1d10 ! This is dummy. True QpGcut_psi is in GWIN0'
      write(ifi,*) ' ------------------------------------------- '
      write(ifi,"(2i4,' ! nbas lmxax (max l for argumentaion)')")
     .nat,lmxax
      write(ifi,*) ' ------------------------------------------- '
      do  isp = 1, nsp
        write(ifi,"(' -- ibas lmxa konf(s) konf(p) konf(d)... '
     .  , ' isp=',2i2)")isp
        do  ib = 1, nbas
          if (lmxa(ib) .gt. -1) then
            write(ifi,"('   ',99i4)")
     .      ipb(ib),lmxa(ib),konf(0:lmxa(ib),ib,isp)
          endif
        enddo
      enddo
      call fclr(' ',ifi)
      end subroutine wlattc


      subroutine mkppovl2(alat,plat,qlat, ng1,ngvec1,ng2,ngvec2,
     inbas, rmax, bas,
     oppovl)

C- < G1 | G2 > matrix where G1 denotes IPW, zero within MT sphere.
c
C     implicit none
      integer ::  nbas, ng1,ng2,nx(3),
     .ig1,ig2, ngvec1(3,ng1),ngvec2(3,ng2),
     .n1x,n2x,n3x,n1m,n2m,n3m
      real(8) :: tripl,rmax(nbas),pi
      real(8) :: plat(3,3),qlat(3,3),
     .alat,tpibaqlat(3,3),voltot, bas(3,nbas)
C     complex(8) :: img =(0d0,1d0)
      complex(8) :: ppovl(ng1,ng2)
      complex(8),allocatable :: ppox(:,:,:)
c-----------------------------------------------------
C     print *,' mkppovl2:'
      pi        = 4d0*datan(1d0)
      voltot    = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      tpibaqlat =  2*pi/alat *qlat
c < G1 | G2 >
      n1x = maxval( ngvec2(1,:)) - minval( ngvec1(1,:))
      n1m = minval( ngvec2(1,:)) - maxval( ngvec1(1,:))
      n2x = maxval( ngvec2(2,:)) - minval( ngvec1(2,:))
      n2m = minval( ngvec2(2,:)) - maxval( ngvec1(2,:))
      n3x = maxval( ngvec2(3,:)) - minval( ngvec1(3,:))
      n3m = minval( ngvec2(3,:)) - maxval( ngvec1(3,:))
c
      allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
      ppox = 1d99
      do ig1  = 1, ng1
        do ig2  = 1, ng2
          nx(1:3) = ngvec2(1:3,ig2) - ngvec1(1:3,ig1) ! G2-G1
          if( ppox(nx(1),nx(2),nx(3)) .eq. 1d99 ) then
            call matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat,
     i      nx(1:3), ! G2 -G1
     o      ppox( nx(1),nx(2),nx(3)))
          endif
        enddo
      enddo
      do ig1 = 1,ng1
        do ig2 = 1,ng2
          nx(1:3) = ngvec2(1:3,ig2) -ngvec1(1:3,ig1) ! G2-G1
          ppovl(ig1,ig2) = ppox( nx(1),nx(2),nx(3) )
        enddo
      enddo
      deallocate(ppox)
      end subroutine mkppovl2


C      subroutine fmain
CC- debugs pvsug1
C      implicit none
C
C      integer nbas
C      parameter (nbas=8)
C      integer ipc(nbas),lmxa(nbas),ipcx(nbas)
C      data ipc/3,1,4,3,2,6,5,1/
C      data lmxa/-1,1,4,-1,2,6,-1,1/
C
C      call pvsug1(nbas,lmxa,ipc,ipcx)
C
C      end
      subroutine pvsug1(nbas,lmxa,ipc,ipcx)

C- Remakes class table, expunging classes with floating orbitals.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   lmxa  :augmentation l-cutoff
Ci   ipc   :class table: site ib belongs to class ipc(ib)
Co Inputs/Outputs
Co   ipcx  :expunged class table: classes with lmxa=-1 are expunged
Co         :and the remaining classes are sequentially renumbered
Co         :preserving the order of the remaining classes.
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Mar 07
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,lmxa(nbas),ipc(nbas),ipcx(nbas)
C ... Local parameters
      integer i,ip,ic,ipskip,nelim
      integer prm(nbas)

      call ivshel(1,nbas,ipc,prm,.true.)
C      do  i = 1, nbas
C        prm(i) = prm(i)+1
C      enddo

C     nelim = number of classes that have been eliminated so far
      nelim = 0
C     ipskip is the number of the last class that was skipped.
C     Multiple occurrences of a skipped class must still only
C     reduce the net number of classes by one.
C     We must avoid incrementing nelim when multiple sites
C     correspond to a skipped class.
      ipskip = 0
C     Loop over sites in order of increasing class index, ip
      do  i = 1, nbas
        ip = prm(i)+1
        ic = ipc(ip)
C       Test whether this class should be purged
        if (lmxa(ip) .lt. 0) then
          if (ipskip .ne. ic) nelim = nelim+1
          ipskip = ic
          ipcx(ip) = -1
        else
          ipcx(ip) = ic - nelim
        endif
      enddo
      end subroutine pvsug1


