      subroutine sugw (jobgw) 
      use m_struc_def,only: s_rv1, s_spec, s_site 
      use m_suham,only: iv_a_oindxo,ham_ldham, ndham=>ham_ndham
      use m_lmfinit, only: 
     &ham_pwmode,ham_pwemin,ham_pwemax,ham_oveps,lrsig=>ham_lsig,stdl,
     &ham_scaledsigma,lat_alat,mxorb,nkaph,nsp,nspc,nl,mxorb,ssite=>v_ssite,sspec=>v_sspec,nbas,lcplxp,n0,nppn,nkap0
      use m_lattic,only: lat_plat, lat_qlat
      use m_supot,only: lat_nabc, lat_gmax
      use m_seneinterp,only: getsenex, senex,dsene
      use m_mksym,only: lat_npgrp,lat_nsgrp
      use m_shortn3,only: shortn3_initialize,shortn3
      use m_mkpot,only: 
     &     smpot=>osmpot, sv_p_osig, sv_p_otau, sv_p_oppi, vconst, vrmt, sv_p_osig, sv_p_otau, sv_p_oppi,
     &     ppn=>ppnl_rv
      use m_mkpot,only: sv_p_osigx,sv_p_otaux,sv_p_oppix,spotx
      use m_MPItk,only: numproc=>nsize,procid,master,master_mpi
      use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_set
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
!! memo: check shortned mechanism.  
!! ----------------------------------------------------------------------
      integer:: i_copy_size,ham_lsig
      integer:: jobgw , lh(n0)
      real(8):: rsml(n0), ehl(n0) 
      logical :: lwvxc,cmdopt0
      integer fopna,fopnx,i,i1,i2,iat,ib,ibr,icore,ierr,ifeigen,!ifi,
     .ifiqg,iflband(2),ifqeigen,ifsyml,igets,igetss,iix,iline, 
     .im1,im2,ipb(nbas),ipqn,ipr,iprint,iq,is,isp,ispc,j,job,k1,
     .k2,k3,konf,konfig(0:n0),l,lchk,ldim,loldpw, ldham(8,2),
     .lgunit,lmaxa,lmxax,lsig,mx,mxint,n1,n2,n3,nat,
     .ncore,ndima,nevl,nev,nglob,ngp,ngp_p, 
     .ngpmx,nline,nlinemax,nlmax,nmx,nn1,nn2,nnn,npgrp,
     .nphimx,npqn,nqbz,nqibz,nqnum,nqnumx,nqtot,nr,nsgrp,iqibz,imx, stdo,
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
      real(8),allocatable::  cphin(:,:,:) !!ovv(:,:),evl_p(:,:), qq1(:,:),qq2(:,:),
      complex(8),allocatable:: ham(:,:),ovl(:,:),evec(:,:),vxc(:,:),
     .ppovl(:,:),phovl(:,:),pwh(:,:),pwz(:,:),pzovl(:,:),
     .testc(:,:),testcd(:),ppovld(:),cphi(:,:,:),cphi_p(:,:,:),
     .geig(:,:,:),geig_p(:,:,:),sene(:,:)
      integer::isize_ham(3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
      integer ndhamx,nspx,nnnx,ifiv
      character strn*120
c#if SIGNN
c      complex(8),allocatable:: sig(:,:)
c#endif
C ... For PW basis
c      integer:: oigv2
c      integer ,allocatable :: kv_iv(:)
      integer pwmode
      double precision pwemin,pwemax,pwgmin,pwgmax,eomin

C ... for reading self-energy
c      integer nqsig
C ... for band plotting
c      real(8),allocatable :: ovvx(:,:,:)
c      real(8) ::ovvpp
c      integer idxdn(n0,nkap0)
c      character lsym(0:n0-1)*1, lorb(3)*1, dig(9)*1, strn4*4
c      data lsym /'s','p','d','f','g','5','6','7','8','9'/
c      data lorb /'p','d','l'/
c      data dig /'1','2','3','4','5','6','7','8','9'/
      integer:: w(1) !! dummy
      real(8):: dnn(3),qlatinv(3,3),qout(3),qtarget(3),qrr(3),axx,bxx,qxxx(3),qxx1(3),qxx2(3)
      integer:: inn(3),iqzz,nqzz,iiiii
      logical:: debug=.false.
      complex(8),allocatable:: evecout(:,:),evecr(:,:)
      real(8),allocatable:: qzz(:,:)
c      integer,allocatable:: igv2(:,:)
      logical:: l_dummy_isanrg,isanrg,oncewrite,newsigmasw,sigmamode !,noshorbz

      integer(4)::   irr,nqi, nmcore !feb2012takao nk1,nk2,nk3, 
c      complex(8):: sfz(nk1,nk2,nk3,ndham,ndham,*) !takao assume ndham =ndimh, OK?

      integer:: nbalance,iqq,nrr,ifiproc,iproc,iadd,ldw
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
      integer:: mode,iwdummy,lso,jx
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),ovlmtoi(:,:),ovliovl(:,:) 
      integer,allocatable:: iantiferro(:)

      integer:: lpdiag=0
      character spid*8

#if MPI | MPIK
      include "mpif.h"
#endif

!! --- Setup ---
      call tcn ('sugw')
      if(iprint()>10) print *,'sugw:'
      sigmamode = mod(lrsig,10) .ne. 0

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
      alat=lat_alat
      n1 = lat_nabc(1)
      n2 = lat_nabc(2)
      n3 = lat_nabc(3)
      plat =lat_plat
      qlat =lat_qlat
      gmax =lat_gmax
      npgrp =lat_npgrp
      nsgrp =lat_nsgrp
      call fftz30(n1,n2,n3,k1,k2,k3)
      stdo = lgunit(1)
      lchk = 1
      lwvxc = .not. cmdopt0('--novxc')

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

c!! Select jobgw
c      if (jobgw .eq. -999) then
c        write(stdo,*) ' lmfgw: Select following job:'
c        write(stdo,*) '  job=-1 : creates files:  GWinput, QPNT, QIBZ, Q0P, QGpsi, QGcou, KPTin1BZ'
c        write(stdo,*) '  job= 0 : init mode; creates files: SYMOPS, LATTC, CLASS, NLAindx, ldima'
c        write(stdo,*) '  job= 1 : GW setup mode; creates: files gwb,gw1,gw2,gwa,vxc,evec,rhoMT.*,normchk'
c        write(stdo,*) '  job=?'
c        read (5,*) jobgw
c      endif
      if(iprint()>10) write(stdo,'(/a,i2)')' gw setup, jobgw', jobgw
      if(jobgw<=-2.or.jobgw>=2) call rx('jobgw number error')
      lsig = 1
c      napw = 0
      i_copy_size=size(ham_ldham)
      call icopy(i_copy_size,ham_ldham,1,ldham,1)
      ldim=ldham(1,1)
      pwmode= ham_pwmode
      pwemin= ham_pwemin
      pwemax= ham_pwemax
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

c$$$!! === Make files SYMOPS, LATTC, CLASS, NLAindx, ldima ===
c$$$      if (jobgw==0 .and. rank0) then
c$$$        call info(30,1,1,' Creating files SYMOPS, LATTC, CLASS, NLAindx, ldima',0,0)
c$$$        allocate(lmxa(nbas),bas(3,nbas)) 
c$$$!! ldima
c$$$        open(newunit= ifildima,file='ldima') 
c$$$        do  ib = 1, nbas
c$$$          lmxa(ib) = sspec(int(ssite(ib)%spec))%lmxa  !int(sspec(ips(ib))%lmxa)
c$$$          if (lmxa(ib) .gt. -1) then
c$$$            call orbl ( ib , 0 , ldim , iv_a_oindxo , norb , ltab , ktab , xx , offl , i1 )
c$$$            write(ifildima,"(3i10)") i1
c$$$          endif
c$$$       enddo
c$$$       close(ifildima)
c$$$!! LATTC
c$$$        lmxax = mxint(nbas,lmxa)
c$$$        allocate(konft(0:lmxax,nbas,nsp))
c$$$        do ib = 1, nbas
c$$$           is = ssite(ib)%spec
c$$$           call dcopy(size(ssite(ib)%pnu), ssite(ib)%pnu,1,pnu,1)
c$$$           call dcopy(size(ssite(ib)%pz),  ssite(ib)%pz, 1,pnz,1)
c$$$           do  isp = 1, nsp
c$$$              do  l  = 0, lmxa(ib)
c$$$                 konft(l,ib,isp) = pnu(l+1,isp)
c$$$                 if( mod(pnz(l+1,isp),10d0) .lt. pnu(l+1,isp) .and. pnz(l+1,isp) .gt. 0) then
c$$$                    konft(l,ib,isp) = mod(pnz(l+1,isp),10d0)
c$$$                 endif    
c$$$              enddo
c$$$           enddo
c$$$        enddo
c$$$        call wlattc(alat,plat,nbas,nat,ipb,lmxax,lmxa,nsp,konft)
c$$$!! NLAindx
c$$$        ifinlaindx = fopnx('NLAindx',2,2,-1)
c$$$        rewind ifinlaindx
c$$$        write(ifinlaindx,'(''----NLAindx start---------------''/I6)') ndima
c$$$        ndima = 0
c$$$!!    This loop order is backwardly compatible with prior versions
c$$$        do  ipqn = 1, 3
c$$$          do  ib = 1, nbas
c$$$            is =  ssite(ib)%spec
c$$$            lmaxa=sspec(is)%lmxa
c$$$            i_copy_size=size(ssite(ib)%pnu)
c$$$            call dcopy  (i_copy_size ,ssite(ib)%pnu,1,pnu,1)
c$$$            i_copy_size=size(ssite(ib)%pnu)
c$$$            call dcopy  (i_copy_size,ssite(ib)%pz, 1,pnz,1)
c$$$            call icopy(sspec(is)%idxdn,sspec(is)%idxdn,1,idxdn,1)
c$$$            if (lmaxa .gt. -1) then
c$$$              do  l = 0, lmaxa
c$$$                npqn = 2
c$$$                if (pnz(l+1,1) .ne. 0) npqn = 3
c$$$                if (ipqn .le. npqn) then
c$$$                  konf = pnu(l+1,1)
c$$$                  if (ipqn .eq. 3) konf = mod(pnz(l+1,1),10d0)
c$$$                  strn4 = dig(konf)//lsym(l)//'_'//lorb(ipqn)
c$$$                  if (idxdn(l+1,1) .eq. 1 .or. idxdn(l+1,2) .eq. 1)
c$$$     .            call chcase(0,1,strn4(2:2))
c$$$                  write(ifinlaindx,'(i6,i3,i4,i6,4x,a)')
c$$$     .            ipqn, l, ipb(ib), ndima, strn4
c$$$                  ndima = ndima + (2*l+1)
c$$$                endif
c$$$              enddo
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$        call fclr(' ',ifinlaindx)
c$$$!! CLASS
c$$$        ificlass = fopnx('CLASS',2,2,-1)
c$$$        rewind ificlass
c$$$        allocate(ipc(nbas),ipcx(nbas))
c$$$        do ib=1,nbas
c$$$           ipc(ib) = ssite(ib)%class !call spackv_array_copy_i8_i('u',ssite(ib)%class,1,ib,ipc)
c$$$        enddo
c$$$        call pvsug1(nbas,lmxa,ipc,ipcx)
c$$$        do  i = 1, nbas
c$$$           if(lmxa(i) .gt. -1) then
c$$$              write(ificlass,'(2I4)') ipb(i), ipcx(i)
c$$$           endif
c$$$        enddo
c$$$        deallocate(ipc,ipcx)
c$$$        call fclr(' ',ificlass)
c$$$        call fexit(0,1,' OK! '//'lmfgw mode=0 generated LATTC CLASS NLAindx ldima',0)
c$$$      endif

!! ==== Read file NLAindx ====
      allocate(nlindx(3,0:lmxax,nat))
      ifinlaindx = fopnx('NLAindx',2,1,-1)
      call ioaindx(3,lmxax,nat,ndima,nlindx,ifinlaindx)
!! ... Read QGpsi !takao feb2012 ---> Not QGcou now
      ifiqg  = fopnx('QGpsi',2,4,-1)
      read(ifiqg ) nqnum, ngpmx ,QpGcut_psi,nqbz,nqi,imx,nqibz
c      endif
!! ==== Write, or read past header information, file gwb ====
      ifigwb = fopna('gwb',-1,4)
      rewind ifigwb
      if(jobgw==1) then
        ef0 = 1d99                !dummy
        ifigwn = fopna('gw1',-1,4)
        rewind ifigwn
        ifigwn = fopna('gw2',-1,4)
        rewind ifigwn
        if(rank0) write(ifigwb) nat,nsp,ndima,ndham,alat,qlat,ef0,nqbz,plat,nqnum,nqi
        allocate(ips(nbas),lmxa(nat),bas(3,nat),iantiferro(nat))
        iat = 0
        do  i = 1, nbas
          lmaxa = int( sspec(int(ssite(i)%spec))%lmxa )
          if (lmaxa .gt. -1) then
            iat = iat + 1
            if (iat .gt. nat) call rx('bug in sugw')
            bas(:,iat)=ssite(i)%pos
            lmxa(iat) = lmaxa
            iantiferro(iat)=ssite(i)%iantiferro
          endif
        enddo
        if(rank0) write(ifigwb) lmxa(1:nat), bas(1:3,1:nat)
        print *,'iantiferro=',iantiferro
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
c          call spacks_copy('u',sspec(is)%name,is,is,spid,is)
          spid=sspec(is)%name
c          print *,' ssssssss spid='//spid//'###'
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
            if(rank0)write(ifigwa) z, nr, a, b, rmt(ib), lmaxa, nsp, ncore,spid
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
        if(rank0) write(ifigwa) iantiferro(1:nat) !iantiferro may2015
        if(rank0) call fclr('gwa',ifigwa)
      else
        call rxi('sugw: bad jobgw',jobgw)
      endif

!! == GW setup loop over k-points ==
      if (lchk .ge. 1 .and. jobgw .eq. 1) then
        ifinormchk = fopna('normchk',-1,0)
        rewind ifinormchk
        write(ifinormchk,849)
  849   format('#     eval          IPW        IPW(diag)    Onsite(tot)   Onsite(phi)      Total')
      endif

!! Write, or read ngpmx
      if(rank0) then
         ifigwb = fopna('gwb',-1,4)
         if (jobgw .eq. 1) then
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
        if (jobgw==1) then
        if (.not. cmdopt0('--novxc')) then
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
         if(ham_scaledsigma/=1d0.and.sigmamode) then
           write(6,*)' Scaled Sigma method: ScaledSigma=',ham_scaledsigma
         endif   
      endif
      iqibz=0

!!  === jobgw divider apr2013takao ===
      nbalance = nqi*nsp/numproc 
      if(mod(nqi*nsp,numproc)/=0) nbalance = nbalance+1
      if(iprint()>10) print *,'nbalance=',nbalance
      if(iprint()>10) write(6,"(' HAM_OVEPS=',d11.2)") ham_oveps
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
      if(pwmode<10) call shortn3_initialize(qlat)
      if(debug)print *,' goto 1001 procid=',procid
      iqq=0
      do 1001 iq = 1, nqtot
        read (ifiqg)  q,  ngp, irr   ! q, and number of G vectors for
        ! eigenfunction expansion at this q
        ! expansion of the Coulomb interaction
        if(irr==1) iqibz=iqibz+1
        lwvxc = lsig>0 .and. jobgw==1 .and. iqibz<=nqibz
        if (cmdopt0('--novxc')) lwvxc = .false.
        if(allocated(ngvecp)) deallocate(ngvecp)
        allocate(ngvecp(3,max(ngp,1)))
        if(ngp==0) then !dummy 
           read (ifiqg)
           ngvecp=-99999
        else
           read (ifiqg)  ngvecp
        endif
        if (jobgw/=1) deallocate(ngvecp) 
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
        qp = q

!!!  ... For this qp, G vectors for PW basis and hamiltonian dimension
        call m_Igv2x_set(qp)    ! get napw and igv2x

        allocate(ham(ndimh,ndimh),ovl(ndimh,ndimh),evec(ndimh,ndimh))
        allocate(vxc(ndimh,ndimh))
        allocate(cphi(ndima,ndimh,nsp),cphin(2,ndimh,nsp))
        do 1002 isp = 1, nsp
          if(iprocq(iqq,isp)/=procid) cycle !cycle when isp for iqq is not for procid 
          lso=0
          if(lso==1) call rx('sugw:lso=1 is not yet')
!! LDA Hamiltonian and overlap matrices for this qp ---
          call hambl(lso,isp,qp, spotx,vconst,sv_p_osigx,sv_p_otaux,sv_p_oppix, vxc,ovl,w) !H(LDA) wihout vxc
          call hambl(lso,isp,qp, smpot,vconst,sv_p_osig, sv_p_otau, sv_p_oppi,  ham,ovl,w) !H(LDA) and O
          vxc = ham - vxc ! vxc(LDA) part 
          if (lwvxc) then
            write(ifiv) ndimh,ldim,qp 
            write(ifiv) q,vxc
          endif
!! LDA + sigma Hamiltonian for this qp ---
          if(sigmamode) then ! See m_bandcal.F
             call getsenex(qp, isp, ndimh,ovl)
             ham(:,:) = ham(:,:) + ham_scaledsigma*senex 
             call dsene() !delete senex
          endif   
!!   --- Branch jobgw = 1 : make cphi, matrix elements ---
          if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
          if (nspc .eq. 2) call rx('diagonalization not ready for nspc=2')
          epsovl = ham_oveps
!! we fined evl(i), z(:,i), i=1,nev !nov2015
          evec=-1d99
          evl(:,isp)=1d99
          call zhev_tk4(ndimh,ham,ovl,ndimh,nev,evl(1,isp),evec,epsovl)
c     if(iq==1)write(stdl,"('fp evl',8f8.4)")(evl(i,isp),i=1,nev)
          write(6,"(' sugw:  kpt isp=',i8,i2,' of ',i8, ' k= ',3f9.5, ' ndimh= ',i5,
     &    ' irank=',i4, ' lwvxc=',l,' nev=',i5)")  iq,isp,nqtot,qp,ndimh,procid,lwvxc,nev
          call prtev(evec,ndimh,evl(1,isp),ndimh,nev)
          if (ndham .gt. nev) evl(1+nev:ndham,isp)=1d20 !this is also given in lmf2gw nov2015. 
          if(mod(iq,10) .ne. 1) call poppr
          if(debug) write(6,"(' sugw:procid iq isp lwvxc= ',3i3,' ',l)")procid, iq,isp,lwvxc 
          if (lwvxc) then
             ifievec = fopna('evec',-1,4)
             write(ifievec) q, evec(1:ndimh,1:ndimh),nev
          endif
          if(magexist) then
             if(isp==1) evl(1:ndimh,isp)=evl(1:ndimh,isp) - vnow/2d0
             if(isp==2) evl(1:ndimh,isp)=evl(1:ndimh,isp) + vnow/2d0
          endif
c$$$  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$  !! wave function rotation test for si 5x5x5.
c$$$  if(.false.) then
c$$$  c       if(iq==8) then
c$$$  ifievv=1012
c$$$  open( ifievv, file='evecx' )
c$$$  ifiqbz=1025
c$$$  open(ifiqbz,file='QBZ')
c$$$  read(ifiqbz,*) nqzz
c$$$  allocate(qzz(3,nqzz))
c$$$  do i=1,nqzz
c$$$  read(ifiqbz,*) qzz(:,i)
c$$$  enddo
c$$$  
c$$$  allocate(evecout(ndimh,ndimh))
c$$$  do iqzz = 1,nqzz
c$$$  qtarget = qzz(:,iqzz)
c$$$  c
c$$$  c                call rotwv(q,qtarget,ndimh,napw,ndimh, plat,qlat,evec,evecout,ierr)
c$$$  call rotwv(q,qtarget,ndimh,napw,ndimh, evec,evecout,ierr)
c$$$  if(ierr/=0) cycle
c$$$  
c$$$  do
c$$$  read(ifievv,*,end=1019) iiiii
c$$$  read(ifievv,*) qrr, ndimh
c$$$  print *,'qrr ndimh=',qrr,ndimh
c$$$  allocate(evecr(ndimh,ndimh))
c$$$  print *,'qrr ndimh xxx=',qrr,ndimh
c$$$  do j= 1,ndimh
c$$$  do i= 1,ndimh
c$$$  read(ifievv,*) axx,bxx
c$$$  evecr(i,j)=dcmplx(axx,bxx)
c$$$  enddo
c$$$  enddo
c$$$  if(sum(abs(qrr-qtarget))<1d-8) exit
c$$$  deallocate(evecr)
c$$$  enddo
c$$$  rewind ifievv
c$$$  c$$$                if(.not.noshorbz()) then
c$$$  c$$$                  call shorbz(qtarget,qxxx,qlat,plat)
c$$$  c$$$                else
c$$$  c$$$                  qxxx=qtarget
c$$$  c$$$                endif
c$$$  qxxx=qtarget
c$$$  write(1013,"(i10)") 11111
c$$$  write(1013,"(3f8.3,i10)") qtarget,ndimh
c$$$  c           write(1013,"(i10,3f8.3,i10)") 11111,qout,ndimh
c$$$  do j=1,ndimh
c$$$  do i=1,ndimh
c$$$  if(abs(evecout(i,j))+abs(evecr(i,j))>1d-4) then
c$$$  write(1013,"(2i4,2d13.5,2x,2d13.5,2x,2d12.4,2x,2d12.4,3f8.3)") 
c$$$  &                i,j,evecout(i,j), evecr(i,j), evecout(i,j)/evecr(i,j),
c$$$  &                abs(evecout(i,j)), abs(evecr(i,j)),qxxx
c$$$  endif
c$$$  enddo
c$$$  write(1013,*)
c$$$  enddo
c$$$  
c$$$  deallocate(evecr)
c$$$  enddo
c$$$  deallocate(evecout)
c$$$  stop 'test end xxxxxxxx'
c$$$  1019         continue
c$$$  stop 'uuuuuuuuuuuuuuu'
c$$$  endif
c$$$  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     !   --- Project wf into augmentation spheres, Kotani conventions ---
c     nlmax = globalvariables%mxorb / globalvariables%nkaph
          nlmax = mxorb / nkaph
          allocate(aus_zv(abs(-nlmax*ndham*3*nsp*nbas)))
          if (-nlmax*ndham*3*nsp*nbas<0) aus_zv(:)=0.0d0
          if(debug) write(6,"(' ppppppp goto makusq')")
          call makusq ( 1 ,  nbas,0, nev,  isp, 1 , qp , evec , aus_zv )
          call gwcphi ( ssite , sspec , isp , nsp , nlmax , ndham , nev 
     .         , nbas , ipb , lmxax , nlindx , ndima , ppn , aus_zv , cphi 
     .         ( 1 , 1 , isp ) , cphin ( 1 , 1 , isp ) )
          if (allocated(aus_zv)) deallocate(aus_zv)
          if(debug) write(6,"(' ppppppp end of gwcphi')")

!     ! We keep note in the followings, but be careful (may contain bugs)...
!     !  --- Overlap of IPWs, PW expansion of eigenfunctions pwz ---
!     !      The IPW basis consisting of PWs with holes knocked out of spheres,
!     !      IPWs must be orthogonalized: IPW -> IPWbar; see PRB 76, 165106.
!     !         |IPWbar_G> = sum_G2 |IPW_G2> O^-1_G2,G1, where O_G1,G2=<IPW_G1|IPW_G2>
!     ! == Definitions ==
C     * ppovl = overlap of IPWs (Generated by pwmat and pwmat2)
C     ppovl_G1,G2 = O_G1,G2 = <IPW_G1 IPW_G2>
C     * pwh = PW expansion of basis function (Generated by pwmat2 only)
C     basis_j> = sum_G2 PWH_G2,j |IPW_G2>
C     * Matrix elements (overlap) of IPW and basis function (Generated by pwmat only)
C     phovl_G1,j = sum_G2 ppovl_G1,G2 pwh'_G2,j  (matrix form PHOVL = O * PWH')
C     Note: phovl is only used as an intermediate construction, old branch
C     * Note: pwh' is expanded to the LMTO cutoff gmax while
C     pwh  is expanded to the GW cutoff QpGcut_psi
C     Thus O^-1 PHOVL will not identically recover PWH.
C     The original branch (loldpw=0) uses effectively PWH'; the new one uses PWH.
C     This is a major distinction between the two (see Remarks)
C     * PW expansion of eigenfunction:
C     |psi_n> = sum_j z_jn |basis_j>
C     * Define pwz_G,n = sum_j PWH_G2,j z_jn  (in matrix form: PWZ = PWH Z)
C     Then
C     |psi_n> = sum_j,G2 z_jn PWH_G2,j |IPW_G2> = sum_G2 PWZ_G2,n |IPW_G2>
C     Overlap of IPW and eigenfunction:
C     PZOVL_G1,n = <IPW(G1) psi_n> = sum_G2 O_G1_G2 PWZ_G2,n
C     PZOVL = O * PWZ (matrix form) <--- old
C     Note: pzovl is only used as an intermediate construction, old branch
          if(debug) print *,'ppppppppppp ngp=',ngp
          if (ngp .gt. 0) then
             allocate(ppovl(ngp,ngp),pwz(ngp,ndimh))
             allocate(phovl(ngp,ndimh))
C     Pass qx to pwmat (or pwmat2):
C     qx = (unshortened) q if 2s digit loldpw = 0
C     qx = (shortened)  qp if 2s digit loldpw = 1
C     Old convention: call pwmat
c     if (mod(loldpw,2) .eq. 0) then

!     !  We have  q+G(igvx; internal in pwmat) = qp + G(igv2)  
!     !  Thus, we have
!     !           igv(internally in pwmat) = igv2 + qlatinv*(qp-q)
!     !            inn = qlatinv*(qp-q)
!!!!  NOTE: from 20Sep2012, we set qp=q (shorbz is in hambl.F. In other words, q supplied to hambl
!!!!  do not need to be short enough).
             if(sum(abs(qp-q))>1d-8) stop 'sugw:qp/=q; qp=p sep2012'
             inn=0
             if(debug) print *,'goto pwmat'
             call pwmat (  ssite , sspec , nbas , ndimh , napw, !slat ,
     .            igv2x, iv_a_oindxo , q , ngp , nlmax , ngvecp , gmax , inn, 
     .            ppovl, phovl )
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
C     --- File output, jobgw=1 ---
          ifigwb = fopna('gwb',-1,4)
          write (ifigwb) q, ndimh, iq
          if(debug) write (6,"('q ndimh=',3f10.5,i10)") q, ndimh
          write (ifigwb) evl(1:ndimh,isp),cphi(:,:,isp)
          write (ifigwb) ngp 
          if(ngp >0 ) write (ifigwb) ngvecp, pwz 
          if(ngp ==0) write (ifigwb) ngvecp
          ifigw1 = fopna('gw1',-1,4)
          write(ifigw1) q, (evl(i,isp),i=1,ndimh)
C     call prmx('e(H+sigma-vxc)',evl,ndimh,ndimh,1)
          
C     ... Overlap checking.   Define:
C     Interstitial part of eigenfunction overlap:
C     <psi_n||psi_n'>
C     = sum_G1,G2 (pwz_G1,n|IPW_G1>)+  (pwz_G2,n'|IPW_G2>)
C     = sum_G1,G2 (pwz_G1,n)+ ppovl_G1,G2 (PWZ_G2,n')
C     = (PWZ)+ O (PWZ) = (PZOVL)+ (PWZ)  (old style)
          if (lchk .ge. 1 .and. ngp .gt. 0) then
             allocate(testc(ndimh,ndimh),testcd(ndimh))
             call zgemm('C','N',ndimh,ndimh,ngp,(1d0,0d0),
     .            pzovl,ngp,pwz,ngp,(0d0,0d0),testc,ndimh)
             deallocate(pzovl)
!     ! call zprm('(PWZ)+^-1 O_i^-1 (PWZ)',2,testc,ndimh,ndimh, ndimh)
             do  i = 1, ndimh
                testcd(i) = sum(dconjg(pwz(:,i))*ppovld*pwz(:,i))
             enddo
             deallocate(ppovld)
C     xx(1) = sum over all augmentation w.f.  cphi+ ovl cphi
C     xx(2) = sum over augmentation phi only.
C     xx(3) = IPW contribution to phi+ phi
C     xx(4) = IPW contribution to phi+ phi, using diagonal part only
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
C     xx(1) = sum(dconjg(cphi(1:ndima,i1,isp))*cphi(1:ndima,i2,isp))
C     xx(2) = sum(dconjg(cphi(1:nchan,i1,isp))*cphi(1:nchan,i2,isp))
                   xx(3) = testc(i1,i2)
                   xx(4) = testcd(i1)
                   if (i1 .eq. i2) then
                      write(ifinormchk,'(f12.5,5f14.6)')
     .                     evl(i1,isp),xx(3),xx(4),xx(1),xx(2),xx(1)+xx(3)
                   endif
                enddo
             enddo
             deallocate(testc,testcd)
             write(ifinormchk,*)
          endif                 !------------- end of lchk=1
          if (ngp .gt. 0) deallocate(ppovl,pwz)
          if (isp .eq. nsp) deallocate(ngvecp) !takao feb2012 ,ngvecc)
!     !   ... Subtract <psi|vxc|psi> from starting evl=e(lda+sigma-vxc)
!     !       evl should contain e(lda+sigma-vxc)
          allocate(testc(ndimh,ndimh))
          call zgemm('C','N',ndimh,ndimh,ndimh,(1d0,0d0),
     .         evec,ndimh,vxc,ndimh,(0d0,0d0),testc,ndimh)
          do i1 = 1, ndimh
             do i2 = 1, ndimh
                evl(i1,isp) = evl(i1,isp) - testc(i1,i2) * evec(i2,i1)
             end do
          end do
          ifigw2 = fopna('gw2',-1,4)
          write(ifigw2) q, (evl(i,isp),i=1,ndimh)
          deallocate(testc)
 1002  continue             ! Loop over spins
        deallocate(ham,ovl,evec,vxc,cphi,cphin)!,igv2)
 1001 continue ! ===== Loop over qp ===============================
      call fclr(' ',ifiqg)
        ifigwb = fopna('gwb',-1,4)
        call fclose(ifigwb)
      deallocate(evl,nlindx)
      call fclr('gw1',fopna('gw1',-1,4))
      call fclr('gw2',fopna('gw2',-1,4))
      call fclr('gwb',fopna('gwb',-1,4))
      if (.not. cmdopt0('--novxc')) then
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
      implicit none
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
      implicit none
C ... Passed parameters
      integer nbas,nat,nsp,lmxax,lmxa(nbas),konf(0:lmxax,nbas,nsp)
      integer ipb(nbas)
      real(8) :: plat(3,3),alat
      integer ifi,ib,isp,fopnx
      open(newunit=ifi,file='LATTC')
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
      close(ifi)
      end subroutine wlattc
c--------------------------------------------------------------------
      subroutine mkppovl2(alat,plat,qlat, ng1,ngvec1,ng2,ngvec2,
     inbas, rmax, bas,
     oppovl)
C- < G1 | G2 > matrix where G1 denotes IPW, zero within MT sphere.
      implicit none
      integer ::  nbas, ng1,ng2,nx(3),
     .ig1,ig2, ngvec1(3,ng1),ngvec2(3,ng2),
     .n1x,n2x,n3x,n1m,n2m,n3m
      real(8) :: tripl,rmax(nbas),pi
      real(8) :: plat(3,3),qlat(3,3),
     .alat,tpibaqlat(3,3),voltot, bas(3,nbas)
C     complex(8) :: img =(0d0,1d0)
      complex(8) :: ppovl(ng1,ng2)
      complex(8),allocatable :: ppox(:,:,:)
      logical:: debug=.false.
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
      if(debug) print *,' mkppovl2: 1 ',n1x,n1m,n2x,n2m,n3x,n3m
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
      if(debug) print *,' mkppovl2: 2 ',n1x,n1m,n2x,n2m,n3x,n3m
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
      implicit none
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


