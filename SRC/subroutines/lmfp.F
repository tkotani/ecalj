      subroutine lmfp(llmfgw) 
      use m_lmfinit,only: irs1,irs2,lhf,bz_nevmx,bz_lmet,bz_w, irs1x10, !,irs3,irs4,irs5,irs11,
     &     ctrl_lfrce,ctrl_mdprm,ctrl_maxit,ctrl_nbas,ctrl_nl,ctrl_nspec,ctrl_nspin,
     &     ctrl_zbak,ctrl_lrel,ctrl_ldos,ctrl_nitmv,ctrl_nvario,ctrl_quit,ctrl_tol,
     &     ctrl_lmet, !,ctrl_lrs,ctrl_lbas
     &     ham_ehf,ham_ehk,ham_seref,ham_eterms,
     &     lat_plat,lat_gam,lat_ldist,lat_alat,lat_dist,lat_plat0,
     &     mix_umix,mix_tolu,mix_b,mix_bv,mix_w,mix_wc,mix_nsave,mix_mmix,
     &     sspec=>v_sspec, ssite=>v_ssite, !v_smix,
     &     sstrnsymg, !ssite2=>v_ssite2,sstrnmix,
     &     symgaf, !v_sbz,v_sctrl,v_sctrl2,v_sham,v_slat,v_slat2,v_spot,v_sarry,v_sarry2,
     &     globalvariables, rv_a_opos 
      use m_supot,only:sv_a_oorhat,zv_a_osmrho,supot
      use m_mksym,only: rv_a_oag,rv_a_osymgr,iv_a_oistab,lat_nsgrp
c      use m_mkqp,only: mkqp
      use m_suham,only:suham
      use m_lattic,only:lattic
      use m_ext,only: sname
!! = Main routine of lmf = (following document is roughly checked at aug2020)
!   lmfp contains two loops
!!    1 big loop do 2000 is for molecular dynamics (relaxiation).
!!    2. innner loop do 1000 is for electronic structure self-consistency.
!!       core of inner loop is bndfp. because of historical reason, Umatrix for LDA+U is managed in lmfp.
!!  All data set to invoke lmfp is stored in m_lmfinit by rdctrl.
!! aug2020. T.kotani removed lshr mode (automatic modification of plat), because
!!      we don't maintain the mode. Probably, we need to re-design it (maybe outside of molecular relaxiation).
!! original classifications are      
C     i   bz_*   :related to the Brillouin Zone
C     i   ctrl_* :               program flow
C     i   ham_*  :               hamiltonian
C     i   pot_*  :               potential
C     i   lat_*  :               lattice
C     i   mix  :                 mixing parameters; see routine umix
!! however, these classification is not so consistent. Some data in m_lmfinit are redundant because
!!      we rewrite old code somehow mechanically.
!!  These are still meaningful, defined in m_struct_def, and set in m_lmfinit
!     i   sspec :struct for species-specific information; see routine uspec
!     i   ssite :struct for site-specific information; see routine usite
!! followings are memo for LDA+U
C     l Local variables
C     l   lmaxu : max l for a U (used for dimensioning)
C     i   lldau(ib)
C     i         : U on site ib with dmat beginning at dmats(*,lldau(ib))
C     l   nlibu : total number of U blocks
Cxxx unused    l   irlxsh: counter for shear relaxations.  irlxsh=0 => new step
C     u Updates
!!  aug2020. T.Kotani remove s_* except s_spec and s_pot. Data are in m_lmfinit.
!!           In future, we need clear protection after data are generated.
C     u   05 Jul 08 Setup for new PW addition to basis
C     u   04 Jul 08 New restart file editor
C     u   20 Jun 06 Repackaged MPI
C     u   21 Mar 06 First cut at shear relaxations
C     u   08 Mar 06 Relaxation restores pos at minimum g when not convgd
C     u   08 Jan 06 can write to restart file rst.#, where # = iter-no.
C     u   09 Nov 05 (wrl) Convert dmat to complex form
C     u   07 Jul 05 rst file version 1.04
C     u   27 Apr 05 LDA+U added (Lambrecht)
C     u   26 Mar 05 Added switch --shorten=no to suppress pos shortening
C     u   23 Feb 05 Bug fix: forces correspondence betw/ pos and site->pos
C     u             after file read of positions.
C     u   11 Jan 05 energy convergence set to ehk when sigma included
C     u   21 Dec 04 Add option to rotate local density on file read
C     u             and to shorten basis vectors after file read
C     u   06 Sep 03 1st cut at automatic optimization of wave function
C     u    9 Jan 03 Undoes lattice shear when writing basis vectors
C     u   21 May 02 Writes restart file after smshft when moving atoms
C     u   15 Feb 02 (ATP) Added MPI parallelization
C     u   08 Jun 01 Revised call to nwit
C     u   15 Feb 01 added density I/O; arguments to bndfp changed.
C     u   17 Jun 00 alpha version.  No relaxations yet.
C     ----------------------------------------------------------------------
c      type(s_sstrn)::sstrn
c      type(s_array)::sarray
c      type(s_bz)::sbz
c      type(s_ctrl)::sctrl
c      type(s_lat)::slat
c      type(s_mix)::smix
c      type(s_ham)::sham
c      type(s_pot)::spot
c      type(s_str)::sstr
c      type(s_move)::smove
c      character  prgnam*8 !sstrn*(*),
c      character*(*)::sstrnmix !,sstrnjobid
c      type(s_spec)::sspec(*)
c      type(s_site)::ssite(*)
      implicit none
      integer ::iwdummy ,idummy
      integer procid,master,mpipid,nproc
      logical mlog
      logical lgors,cmdopt,bittst
      integer fopna,fopng,i,ifi,igets,iors,ipr,irs1x,isw,ix(5),j, !,iscr
     .     k,lcgf,leks,lgunit,lpnu,lrel,lrout,nbas,nat,
     .     nbaspp,nevmx,nglob,nit1,nl,nsp,nspec,numq, !npadl,npadr,lsx,
     .     stdo,pdim,lsc
      integer:: o
      real(8),allocatable  :: pos_move(:,:)
      real(8) ,allocatable :: frc_rv(:)
      real(8) ,allocatable :: ftot_rv(:)
      real(8) ,allocatable :: pos2_rv(:,:)
      real(8) ,allocatable :: wk_rv(:)
      real(8) ,allocatable :: p_rv(:)
      double precision plat(3,3),qlat(3,3),qbg,xv(10),fptol,umix
      character strn*120, fileid*68, alabl*8, flg*3

C     For mixing.  Default parameters dmxp:
C     1(I) mixing scheme; 2 beta; 3 wc; 4,5 wgt; 6(I)nsave 7(I) mmix;
C     8(I) nkill; 9 betv; 10 rmscst.  11-20 are outputs:
C     11 rmsdel 12 rms2 13 nmix 14 actual broy 15 actual beta 16-24 tj
C     25 1 if wt ne 0, 10 if wa ne 0, 11 if all nonzero
C     27..29: hold parms for static parms block regular mixing
C     30..32: hold parms for static parms block Euler angle mixing
C     33 : Lindhard screening parameter
      double precision dmxp(33)
C     ... for iterations
      logical lbin,a2bin !lhf,
      integer maxit,iter
      double precision seref,etot(2),amom,qdiff,qtol,etol,alat
c      equivalence (qdiff,dmxp(11))
C     ... for relaxation
      logical ::xyzfrz(3),ltmp !,lshr=.false.
      integer:: icom , natrlx , nvrelx , ltb , itrlx , nm , irlxsh 
     .     , nitrlx , bitor , ng
      integer ,allocatable :: indrx_iv(:)
      real(8) ,allocatable :: w_rv(:)
c     ki, for rlse and allocate oindrx
      integer,allocatable:: iv_tmp(:)
      double precision mdprm(6),gam(4),gam1,bstim,rhosig,pletot(6,2),
     .     plat0(3,3),dist0(9),dist(9)
      parameter (nm=3)
C     ... for LDA+U
      integer nlibu,lmaxu
c     integer:: odmatu , ovorbdmat , odmato
      complex(8),allocatable::vorbdmat(:),dmatu(:),dmato(:)
      integer ,allocatable :: lldau_iv(:)
      double precision tolu
      data irlxsh /0/ dist0 /9*0d0/
      logical:: l_dummy_isanrg,isanrg
      real(8):: pnu(20),bz_wx
      logical :: llmfgw,nowriteden=.false., nexist
      integer:: isize_iv_tmp,ibas,unlink,ifipos,ifile_handle, ifivorbdmat,nsize

#if (MPI |MPIK)
      integer:: ierr
      include "mpif.h"
#endif
      call tcn('lmfp')
      etot(1) = 0
!! MPI-specific
      nproc  = mpipid(0)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
!!  -------------------------- Total energy mode -->moved to m_lmfinit ---
c$$$      if (cmdopt('--etot',6,0,strn)) then
c$$$!!  No forces or dynamics
c$$$        mdprm=0
c$$$        ctrl_lfrce=0
c$$$        ctrl_mdprm=0
c$$$        nowriteden=.true.
c$$$        ctrl_maxit=1
c$$$      endif
!!  initialization -------------------
      call getpr(ipr)
      nbas=ctrl_nbas
      nl=ctrl_nl
      nspec=ctrl_nspec
      nsp=ctrl_nspin
      nat = globalvariables%nat
      qbg = ctrl_zbak(1)
      maxit = int(ctrl_maxit)
      lrel = isw ( int(ctrl_lrel) .ne.0 )
      call setcc(lrel) !lrel=T means scaler relativistiv c=274.074d0 in a.u.
      if (lhf) maxit = 1
      nbaspp = nbas      
      stdo   = lgunit(1)
      irs1x = irs1
      plat = lat_plat
      if (ipr .ge. 30) call praugm(sspec,0) !Printout properties of species
      
!! ... --rs=3 => always read from atom file !2020aug.
!! Sep2020 takao: --rs=3 mode is removed. (--rs=3 meand fixed density Harris-foukner MD).
c      if ( iand(3,int(ctrl_lrs)) .eq. 3 ) irs1x = 0 
C     ... Setup for no screening transformation
!! sep2020 comment out below      
c$$$C     Shorten site positions
c$$$      if (.not. cmdopt('--shorten=no',12,0,strn)) then
c$$$         allocate(pos2_rv(3,nbas))
c$$$         do ibas=1,nbas
c$$$            pos2_rv(:,ibas) = ssite(ibas)%pos
c$$$         enddo
c$$$         ix(1) = 2
c$$$         ix(2) = 2
c$$$         ix(3) = 2
c$$$         call info0(50,1,0,' lmfp : shortening basis vectors ... ')
c$$$         call shorps ( nbas , plat , ix , pos2_rv , rv_a_opos )
c$$$         do ibas=1,nbaspp
c$$$            ssite(ibas)%pos = rv_a_opos(:,ibas)
c$$$         enddo
c$$$         if (allocated(pos2_rv)) deallocate(pos2_rv)
c$$$      endif

!! ... Setup for charge mixing
      dmxp=0d0
      dmxp(2)=mix_b
      dmxp(9)=mix_bv
      dmxp(4:5)= mix_w(1:2)
      dmxp(3)= mix_wc
      dmxp ( 6 ) = int(mix_nsave)
      dmxp ( 7 ) = int(mix_mmix)
      call parms0(0,0,0d0,0)
C     ... Allocate memory for forces
c      lfrce = int(ctrl_lfrce)
      if (ctrl_lfrce .ne. 0) then
        numq = 1
        if ( int(bz_lmet) .eq. 4 ) numq = 3
        allocate(frc_rv(3*nbas*numq))
      endif
C     ... Relaxation setup
      itrlx = 1
C     Initial shear was already folded into plat
      gam=lat_gam
      gam1 = gam(4)
      gam(4) = 1
      lat_gam=gam
      lat_ldist=0
      nitrlx=ctrl_nitmv
      mdprm=ctrl_mdprm
c      lshr = nint(mdprm(1)) .gt. 100
      if (nint(mdprm(1)) .eq. 0) nitrlx = 0
      if (nint(mdprm(1)) .gt. 0 .and. nint(mdprm(1)) .lt. 4) then
        call rx('lmf not set up for MD yet')
      endif
      if (nitrlx > 0) then
        ctrl_mdprm=mdprm
        allocate(indrx_iv(6*nbas))
C     Next lines in case lattice relaxation
c$$$        if (lshr) then
c$$$          if (abs(gam(4)-1) .gt. 1d-10) call rx('lmfp: '//
c$$$     .           'use of SHEAR= incompatible w/ lattice relaxation')
c$$$          plat0=lat_plat0
c$$$        endif
        call rlxstp (  ssite , natrlx , nvrelx , indrx_iv ,  xyzfrz , pdim )
        icom = 0
        if (nvrelx .ne. 0) then
C     ki The new indrx_iv uses the content of the old indrx_iv.
          isize_iv_tmp=size(indrx_iv)
          allocate( iv_tmp(isize_iv_tmp) )
          iv_tmp(:)= indrx_iv(:)
          deallocate(indrx_iv)
          allocate(indrx_iv(2*natrlx))
          isize_iv_tmp=min(isize_iv_tmp,2*natrlx)
          indrx_iv(:isize_iv_tmp)=iv_tmp(:isize_iv_tmp)
          deallocate(iv_tmp)
          allocate(w_rv(nvrelx*nvrelx))
          allocate(p_rv(pdim))
        endif
        alat = lat_alat
        if(procid==master) then
           open(newunit=ifipos,file='AtomPos.'//trim(sname),form='unformatted',status='new')
           write(ifipos) nbas
           write(ifipos) itrlx,rv_a_opos
        endif   
c$$$        if (procid .eq. master) then
c$$$          ifi = fopna('bsmv',-1,0)
c$$$          allocate(pos2_rv(3,nbas))
c$$$          j = 1
c$$$          call ivset(ix,1,3,j)
c$$$          call shorps ( nbas , plat , ix , rv_a_opos , pos2_rv )
c$$$          call iobsm0 ( 0 , bstim , 0d0 , 0d0 , nbas , alat , pos2_rv, ifi )
c$$$          if (allocated(pos2_rv)) deallocate(pos2_rv)
c$$$        endif
      endif
      
!!=== Re-entry for lshr (share mode); share mode is currently commented out.
!!      Probably share mode is going to be moved to outside of lmfp =======
c    4 continue !
      write(stdo,*)' lmfp : potential and basis setup ... '
      call supot(0) 
      call suham(sspec,ssite)
      if ( ctrl_quit .eq. 8 ) call rx0('quit = ham')
      allocate(ftot_rv(3*nbas))
      nevmx = bz_nevmx !Maximum number of eigenvalues
      lrout = 1  ! Whether to evaluate output density and/or KS energy
      leks = 1
      j = 6
      if (nevmx .eq. -1) then
        lrout = 0
        leks = 0
        bz_nevmx=nevmx
      endif
      lpnu = 1 !Whether to float pnu's
      if (lrout == 0 .and. ctrl_lfrce /= 0) then
        write(stdo,"(a)") 'lmfp (fatal): output density required when forces sought.\n'//
     &              '      To make output density turn off HF=t and/or NEVMX<0'
        call rx('incompatible input. see the end of console output')
      endif
      if (lrout == 0 .and. cmdopt('--etot',6,0,strn)) then
        write(stdo,"(a)") 'lmfp (fatal): output density required with --etot switch.\n'//
     &              '      To make output density turn off HF=t and/or NEVMX<0'
        call rx('incompatible input. see the end of console output')
      endif
      if (lrout == 0 .and. maxit > 1) then
         write(stdo,"(a,i5,a)") 'lmfp (warning): ',maxit,
     &    ' iterations sought but no output rho ... do 1 iteration'
        maxit = 1
      endif
!!... LDA+U initialization
      allocate(lldau_iv(nbas))
      lldau_iv(:)=0 
      call suldau ( nbas , sspec , ssite , nlibu , lmaxu , lldau_iv )
      if (nlibu .gt. 0) then !! ... return nlibu > 0 if any U blocks
        i = nsp*nlibu*(lmaxu*2+1)**2
      else
        i = 1
      endif
      allocate(vorbdmat(i),dmatu(i),dmato(i))
      vorbdmat = 0d0
      dmatu = 0d0
      dmato = 0d0
      if (nlibu > 0) then ! need group info to symmetrize site density matrix
        ng=lat_nsgrp
        umix=mix_umix
        tolu=mix_tolu
        if (umix .eq. 0) umix = 1
        call sudmtu ( nbas , nsp , nlibu , lmaxu , ssite , sspec , 0
     .        , lldau_iv , ng , rv_a_osymgr , iv_a_oistab , dmatu , vorbdmat)
        if ( int(ctrl_quit) .eq. 16 ) call rx0('quit = dmat')
      endif
     
!! === do loop for molecular dynamics =======
      if(nitrlx>=1) allocate(pos_move(3,nbas))
      do 2000 itrlx   = 1,max(1,nitrlx)
!! ===   do loop for electronic structure ===
         do 1000 iter = 1,max(1,maxit)
c     5 continue
C     --- Read restart file or overlap free atom densities ---
C     irs(1) tells what to read and whether to invoke smshft.
C     0+1's bits irs(1)     action
C     0           read from atom file
C     1             read from binary rst file
C     2             read from ascii rsta file
C     3             read nothing (data already input)
C     4s' bit of irs(1) -> invoke smshft after file read.
C     8s' bit of irs(1) -> rotate local density after file read
c 10         continue
            if (irs1x == 0) then !read overlap free-atom densities
               call rdovfa ( nbas, nspec, ssite, sspec,qbg, sv_a_oorhat ) !initial potential
               nit1 = 0
            elseif ( 1<= mod(irs1x,4) .and. mod(irs1x,4)<=2 ) then 
               lbin = irs1x==1 !.not. bittst(irs1x,2) !lbin=T:rst --rs=1, or lbin=F:rsta --rx=2,
               k = -1
               if (procid .eq. master) then
                  ifi=ifile_handle() 
                  if(lbin)       open(ifi,file='rst.'//trim(sname),form='unformatted')
                  if(.not. lbin) open(ifi,file='rsta.'//trim(sname)) 
               endif
               call mpibc1(ifi,1,2,mlog,'lmfp','ifi')
               bz_wx = bz_w
               k = iors ( ssite , sspec ,  fileid 
     .              , nbas , nat , nspec , sv_a_oorhat , iwdummy , nit1 , lbin , ifi )
               bz_w = bz_wx
               if (procid .eq. master) close(ifi) !call fclose(ifi)
               call mpibc1(k,1,2,mlog,'lmfp','k')
               do ibas=1,nbaspp 
                  rv_a_opos(:,ibas) = ssite(ibas)%pos
               enddo
               if (k .lt. 0) then
                  irs1x = 0
                  call rdovfa ( nbas, nspec, ssite, sspec,qbg, sv_a_oorhat ) !initial potential
                  nit1 = 0
c                  goto 10
               elseif (irs1x10) then  
                  call smshft(ssite, sspec, sv_a_oorhat, zv_a_osmrho, 1)
               endif
            endif               !!if not, read nothing because of previous iteration
!! ... Write positions after file read, and repack
            if (ipr .ge. 50) then
               write(stdo,357) 'Basis, after reading restart file'
 357           format(/1x,a/' site spec',8x,'pos (Cartesian coordinates)',9x,
     .              'pos (multiples of plat)')
               call dinv33(plat,1,qlat,xv)
               do  i = 1, nbas
                  j=ssite(i)%spec
                  xv(1:3)=ssite(i)%pos
                  alabl=sspec(j)%name
                  call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
                  write(stdo,345) i, alabl, (xv(j),j=1,3), (xv(3+j),j=1,3)
 345              format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
               enddo
            endif
c$$$  C     --- Optionally re-shorten basis vectors ---
c$$$  if (cmdopt('--shorps',8,0,strn)) then
c$$$  allocate(pos2_rv(3,nbas))
c$$$  do ibas=1,nbas
c$$$  pos2_rv(:,ibas ) = ssite(ibas)%pos 
c$$$  enddo
c$$$  ix(1) = 2
c$$$  ix(2) = 2
c$$$  ix(3) = 2
c$$$  call shorps ( - nbas , plat , ix , pos2_rv , rv_a_opos )
c$$$  call info0(20,1,-1,' lmfp  : write shortening vectors to file shorps ...')
c$$$  call iopos ( .true. , 1 , 'shorps' , nbas , rv_a_opos )
c$$$  call shorps ( nbas , plat , ix , pos2_rv , rv_a_opos )
c$$$  do ibas=1,nbaspp
c$$$  ssite(ibas)%pos= rv_a_opos(3*(ibas-1)+1: 3*(ibas-1)+3)
c$$$  enddo
c$$$  if (allocated(pos2_rv)) deallocate(pos2_rv)
c$$$  endif


C     Hang on to previous site density matrix for this iteration
            if (nlibu .gt. 0) then
               dmato=dmatu
               dmatu=0d0
            endif
!!
            inquire(file='vorbdmat.dat',exist=nexist,size=nsize)
            if(iter==1.and.nexist.and.nsize>0) then
               open(file='vorbdmat.dat',newunit=ifivorbdmat,form='unformatted')
               read(ifivorbdmat) vorbdmat
               close(ifivorbdmat)
            endif

C     --- Make and diagonalize hamiltonian, make new charge density ---
            if (maxit .eq. 0) call info0(20,1,0,' lmfp  : zero iterations sought ... no band pass')
            call bndfp ( nlibu , lmaxu , lldau_iv , ssite , sspec ,    
     .           leks , lrout , ctrl_lfrce , lpnu , dmxp , iter , maxit  , !evl_rv, 
     .           frc_rv , dmatu, vorbdmat, llmfgw )
C     ... check convergence of dmatu and update it and vorbdmat if necessary
            if (nlibu .gt. 0 .and. maxit .gt. 0 .and. lrout .gt. 0) then
               call chkdmu ( nbas , nsp , nlibu , lmaxu , ssite , sspec !, sham 
     .              , 0 , dmatu , dmato, vorbdmat, tolu , umix
     .              , lldau_iv , ng , rv_a_osymgr , iv_a_oistab )
               if(procid==master) then
                 open(file='vorbdmat.dat',newunit=ifivorbdmat,form='unformatted')
                 write(ifivorbdmat) vorbdmat
                 close(ifivorbdmat)
               endif
            endif

!! Write smoothed charge density for contour plotting ---
!!   --density is now moved to locpot.F(rho1mt and rho2mt) and mkpot.F(smooth part).
C     --- Write restart file (skip if --quit=band) ---
            if (procid .eq. master) then
               if ( ctrl_quit .ne. 4 ) then
C     Suppress saving rst file in the middle of a shear (irlxsh > 0)
                  if (irs2.gt.0 .and. (lrout.gt.0 .or. maxit .eq. 0) .and. irlxsh .eq. 0) then
                     lbin = irs2 .ne. 2
c$$$                     if (lbin) fileid = 'rst'
c$$$                     if (.not. lbin) fileid = 'rsta'
c$$$                     if (irs2 .eq. 3) then
c$$$                        call word(fileid,1,i,j)
c$$$                        j = j+1
c$$$                        fileid(j:j) = '.'
c$$$                        call bin2a(' ',0,0,iter,2,0,len(fileid),fileid,j)
c$$$                        if (lbin) ifi = fopng(fileid,-1,8+4)
c$$$                        if (.not. lbin) ifi = fopng(fileid,-1,8)
c$$$                        call info0(10,1,-1,' lmfp:  writing to restart file '//fileid)
c$$$                     else
c$$$                        if (lbin) ifi = fopna(fileid,-1,4)
c$$$                        if (.not. lbin) ifi = fopna(fileid,-1,0)
c$$$                     endif
                  ifi=ifile_handle() 
                  if(lbin)       open(ifi,file='rst.'//trim(sname),form='unformatted')
                  if(.not. lbin) open(ifi,file='rsta.'//trim(sname)) 

                  fileid = 'lmfp:  ' 
                  k = iors ( ssite , sspec ,   fileid !sbz ,, sctrl slat , spot , 1  ,
     .                    , nbas , nat , nspec , sv_a_oorhat , iwdummy , iter , lbin , 
     .                    - ifi )
                  close(ifi) !call fclose(ifi)
                  endif
               endif
            endif
            if (cmdopt('--window=',9,0,strn)) call rx0('lmf : early exit (--window option)')
C     --- Add to save file; decide on next iteration ---
            if (maxit <= 0) goto 9998
            etot(1)=ham_ehf
            etot(2)=ham_ehk
            seref = ham_seref   !   ... Subtract reference energy
            etot(1) = etot(1) - seref
            if (etot(2) .ne. 0) etot(2) = etot(2) - seref
            amom = ham_eterms(15)
C     The desired tolerances in q,e
            qtol = ctrl_tol(1)
            etol = ctrl_tol(3)
#if MPI|MPIK
            flush(6)
            call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
            if (procid .eq. master) then
               rhosig = ham_eterms(19)
               i = 0
               if (rhosig .ne. -99 .and. rhosig .ne. 0) i = 10
C     kino nwit() Outputs
C     kino   lsc   :0 self-consistency achieved (diffe<=etol, qdiff=dmxp(11)<=qtol)
C     kino         :1 if not self-consistent, but encountered max. no. iter.
C     kino         :2 Harris energy from overlap of free atoms (iter=1 and lhf=t)
C     kino         :3 otherwise
               call nwit ( int(ctrl_nvario) , iter , maxit , lhf.or.irs1x==0.and.iter==1,
     &              leks+i , etol , qtol , dmxp(11) , 'cxhi' , amom , etot , lsc ) !dmxp(11)=qdiff 
            endif
            call mpibc1(lsc,1,2,mlog,'lmfp','lsc')
            if (lsc .eq. 2 .and. .not. lhf .and. maxit .gt. 1) lsc = 3
            if (lsc .eq. 1 .and. lrout .gt. 0  .or. lsc .eq. 3) then
               if (iter .ge. maxit) lsc = 1
               if (iter .lt. maxit) lsc = 3
            endif
            if ( int(ctrl_quit) .eq. 4 ) call rx0 ( 'lmf : exit (--quit=band)' )
C     Continue iterations toward self-consistency
            if (lsc .gt. 2) then
               irs1x=3
               cycle            !goto 5
            endif
            exit
 1000    continue ! ---------------- SCF (iteration) loop end ----
         if (nitrlx==0) exit
!! --- Molecular dynamics -----------------------------
         do ibas=1,nbas
            ssite(ibas)%pos0 = ssite(ibas)%pos
         enddo   
         do ibas=1,nbas
            pos_move(:,ibas) = ssite(ibas)%pos
         enddo
         mdprm=ctrl_mdprm
c$$$  if (lshr) then !shear mode
c$$$  call grdepl ( nvrelx , indrx_iv , 0.01d0 , etot(1) , irlxsh , 
c$$$  .           pletot , dist )
c$$$  if (irlxsh .ne. 0) then
c$$$  call grdep2 ( 1 , nvrelx , indrx_iv , dist0 , dist )
c$$$  call newshear()
c$$$  if ( iand(3,int(ctrl_lrs)) .eq. 3 ) then !!! Decide on what density to use
c$$$  irs(1) = 0
c$$$  else             !  Else, use file density
c$$$  irs(1) = iand(7,int(ctrl_lrs))
c$$$  endif
c$$$  goto 4
c$$$  else
c$$$  call relax ( ssite , sspec , itrlx , indrx_iv !sctrl ,
c$$$  .              , natrlx , nvrelx , pletot , p_rv , w_rv , 0 , 0d0 , dist0 
c$$$  .              , icom )
c$$$  call dpzero(dist,6)
c$$$  call grdep2 ( 1 , nvrelx , indrx_iv , dist0 , dist )
c$$$  dist(7) = 1
c$$$  endif
c$$$  else !without shear
         call relax ( ssite , sspec , itrlx , indrx_iv !sctrl ,
     .        , natrlx , nvrelx , frc_rv , p_rv , w_rv , 0 , 0d0 , pos_move 
     .        , icom )
c$$$  endif

!! we comment our next procedure at sep2020  for simplification
c$$$C     Restore lattice symmetry to machine precision
c$$$         if (cmdopt('--fixpos',8,0,strn)) then
c$$$            ng=lat_nsgrp
c$$$            j = 8+1
c$$$            if (strn(9:13) .eq. ':tol=') j = 13
c$$$            if (strn(9:9) .ne. ':' .or. .not. a2bin(strn,fptol,4,0,' ',j,len(strn))) fptol = 1d-5
c$$$            call fixpos ( pos_move , nbas , fptol , ng , plat , rv_a_osymgr, rv_a_oag , iv_a_oistab )
c$$$         endif
         
c$$$  C     Write updated positions to bsmv file
c$$$  if (procid .eq. master) then ! .and. .not. lshr) then
c$$$  ifi = fopna('bsmv',-1,0)
c$$$  call poseof(ifi)
c$$$  bstim = bstim+1
c$$$  allocate(pos2_rv(3,nbas))
c$$$  j = 1
c$$$  call ivset(ix,1,3,j)
c$$$  call shorps ( nbas , plat , ix , pos_move , pos2_rv )
c$$$  call iobsmv ( 0 , bstim , 0d0 , 0d0 , nbas , alat , pos2_rv 
c$$$  .           , - ifi )
c$$$  if (allocated(pos2_rv)) deallocate(pos2_rv)
c$$$  call fclose(ifi)
c$$$  endif
!!  Updated positions in site structure !WARN! t.kotani think this is confusing because
!!  CAUTION: 'positions written in ctrl' and 'positions written in rst' can be different.
         if(procid==master) write(ifipos) itrlx,pos_move
         do ibas=1,nbas
            ssite(ibas)%pos = pos_move(:,ibas)
         enddo
C     ... Exit when relaxation converged or maximum number of iterations
         if (icom .eq. 1) then
            if (procid .eq. master) then
               k = int(ctrl_nvario)
               flg = 'C67'
               call nwitsv(1+2,k,flg,nsp,amom,etot)
            endif
            if ( procid==master) then
               call tcx('lmfp')
               call fexit(0,111,
     .              ' LMFP: relaxation converged after %i iteration(s)',itrlx)
            else
               call tcx('lmfp')
               call fexit(0,111,' ',0)
            endif
c     else
c     call query('proceed with next relaxation step',-1,0)
         endif

C     ... Restore minimum gradient positions if this is last step
C     if (itrlx .eq. nitrlx .and. icom .eq. -1) then
         if (itrlx .eq. nitrlx) then
c     if (.not. lshr) then
            call info0(20,1,0,' lmfp: restore positions for minimum g')
C     call prmx('initial positions',w(opos),3,3,nbas)
            call prelx1 ( 1 , nm , .false. , natrlx , nvrelx , indrx_iv , p_rv !lshr=F
     .           , pos_move )
C     call prmx('minimum-g positions',w(opos),3,3,nbas)
c$$$  else
c$$$  call prelx1 ( 1 , nm , lshr , natrlx , nvrelx , indrx_iv , 
c$$$  .              p_rv , dist0 )
c$$$  call dpzero(dist,6)
c$$$  call grdep2 ( 1 , nvrelx , indrx_iv , dist0 , dist )
c$$$  
c$$$  call info2(20,0,0,
c$$$  .              ' lmfp : strain of minimum gradient:'//
c$$$  .              '%N   PDEF=%6;8,4D'//
c$$$  .              '%N STRAIN=%6;8,4D',
c$$$  .              dist0,dist)
c$$$  endif

C     Repack updated positions in site structure
            do ibas=1,nbas
               ssite(ibas)%pos = pos_move(:,ibas)
            enddo
         endif
!     ! New density after atom shifts 
!     ! If explicitly told to read from atom files after atom movmment
c     if ( iand(3,int(ctrl_lrs)) .eq. 3 ) then
c     irs1x = 0
c     !     ! Else, use self-consistent
c     else !if (.not. lshr) then
         irs1x=3
         call smshft(ssite,sspec,sv_a_oorhat,zv_a_osmrho,ctrl_lfrce) !sctrl,sham,slat,
c     endif
!     ! Write restart file (to include new positions)
         if (procid .eq. master) then ! .and. .not. lshr) then
            ifi = ifile_handle()
            open(ifi,file='rst.'//trim(sname),form='unformatted') !fopna('rst',-1,4)
            fileid = 'lmfp:  '  !// trim(sstrnjobid) !sstrn(i:j)
            k = iors ( ssite , sspec ,  fileid !sbz ,, sctrl slat , spot , 1 , 
     .           , nbas , nat , nspec , sv_a_oorhat , iwdummy , iter , .true. 
     .           , - ifi )
            close(ifi) !call fclose(ifi)
         endif
         if (procid .eq. master) then
            write(stdo,*)' Delete mixing and band weights files ...'
            ifi = unlink('mixm.'//trim(sname))
            ifi = unlink('wkp.'//trim(sname))
         endif
         call parms0(0,0,0d0,0) !   reset mixing block
!     ! Exit when maximum number of iterations encountered
         if (itrlx .eq. nitrlx) then
            if (procid .eq. master) then
               call tcx('lmfp')
               call fexit(1,111,' LMFP: relaxation incomplete after %i iteration(s)',nitrlx)
            else
               call tcx('lmfp')
               call fexit(1,111,' ',0)
            endif
         endif
c     itrlx = itrlx+1
c$$$  if (lshr) then
c$$$  call newshear()
c$$$  !! Decide on what density to use
c$$$  if ( iand(3,int(ctrl_lrs)) .eq. 3 ) then
c$$$  irs(1) = 0
c$$$  else                !  Else, use file density
c$$$  irs(1) = iand(7,int(ctrl_lrs))
c$$$  endif
c$$$  goto 4
c$$$  else
c     goto 5
c     endif
c     cycle
c     c      endif                     ! if (nitrlx .gt. 0 .and. lsc .le. 2) then
 2000 continue
      
!! Write positions to file
 9998 continue
      if (cmdopt('--wpos=',7,0,strn) .or. cmdopt('--wpos:mode1:',13,0,strn)) then
         allocate(pos_move(3,nbas))
         do ibas=1,nbas
            pos_move(:,ibas) = ssite(ibas)%pos
         enddo
         gam=lat_gam
         allocate(wk_rv(3*nbas))
         call rdistn (pos_move , wk_rv , nbas , gam (1) , gam(2) , gam ( 3 ) , 1 / gam1 )
         if (allocated(wk_rv)) deallocate(wk_rv)
         if (allocated(pos_move)) deallocate(pos_move)
      endif 
      if(procid==master.and.nitrlx>0) close(ifipos)
      call tcx('lmfp')
      if (allocated(lldau_iv)) deallocate(lldau_iv)
      if (allocated(ftot_rv)) deallocate(ftot_rv)
      if (allocated(p_rv)) deallocate(p_rv)
      if (allocated(w_rv)) deallocate(w_rv)
      if (allocated(indrx_iv)) deallocate(indrx_iv)
      if (allocated(frc_rv)) deallocate(frc_rv)
      return
      end subroutine lmfp

c$$$      contains
c$$$      subroutine newshear()
c$$$! Store lat_plat,lat_ldist,lat_dist, lattic, mkqp
c$$$! not tested well yet.      
c$$$      if (procid .eq. master) then
c$$$        call info0(20,0,0,' Delete mixing and band weights files ...')
c$$$        ifi = fopna('mixm',-1,4)
c$$$        call dfclos(ifi)
c$$$        ifi = fopna('wkp',-1,4)
c$$$        call dfclos(ifi)
c$$$      endif
c$$$!!  Restore plat, pos to their undistorted state:
c$$$!!  undo original transformation = P P_0^-1
c$$$      allocate(pos2_rv(3,nbas))
c$$$      call dinv33(plat,0,xv,xv(10))
c$$$      call dgemm ( 'N' , 'N' , 3 , nbas , 3 , 1d0 , xv , 3 , pos_move 
c$$$     .     , 3 , 0d0 , pos2_rv , 3 )
c$$$      call dgemm ( 'N' , 'N' , 3 , nbas , 3 , 1d0 , plat0 , 3 , pos2_rv
c$$$     .     , 3 , 0d0 , rv_a_opos, 3 )
c$$$      do ibas=1,nbas
c$$$        ssite(ibas)%pos  = pos_move(3*(ibas-1)+1: 3*(ibas-1)+3)
c$$$      enddo
c$$$      call cppos(1,nbas,ssite)
c$$$!! New shear
c$$$      lat_plat=plat0
c$$$      lat_ldist=3
c$$$      lat_dist=dist
c$$$!! A little memory leakage rel to 1st pass, but not so serious
c$$$      if (allocated(pos2_rv)) deallocate(pos2_rv)
c$$$      if (allocated(lldau_iv)) deallocate(lldau_iv)
c$$$      if (allocated(ftot_rv)) deallocate(ftot_rv)
c$$$      deallocate(vorbdmat,dmato,dmatu)
c$$$      call lattic() !nbas) !slat,
c$$$      plat=lat_plat
c$$$!! Remake qp (t.kotani wonders this mkqp is needed or not?
c$$$      ltmp = iand(1,int(ctrl_lmet)) .ne.0 .or. iand(4+2+1,int(ctrl_ldos)) .ne.0
c$$$      call mkqp(ltmp,.false.,1,-2) 
c$$$!! Write restart file (to include new positions)
c$$$      if (procid .eq. master .and. irlxsh .eq. 0) then
c$$$        ifi = fopna('rst',-1,4)
c$$$        fileid = 'lmfp:  ' !// trim(sstrnjobid) !sstrn(i:j)
c$$$        k = iors (ssite , sspec ,   fileid !sbz ,, sctrl slat ,spot ,  1  , 
c$$$     .        , nbas , nat , nspec , sv_a_oorhat , iwdummy , iter , .true. 
c$$$     .        , - ifi )
c$$$        call fclose(ifi)
c$$$      endif
c$$$      end subroutine newshear
c$$$      end subroutine lmfp



c$$$      subroutine cppos(ib1,ib2,ssite)
c$$$      use m_struc_def           !Cgetarg
c$$$C     - Copy site positions to p0 for a range of sites
c$$$C     implicit none
c$$$      integer ib1,ib2
c$$$      real(8):: pos(3)
c$$$      type(s_site)::ssite(*)
c$$$      integer ib
c$$$      do  ib = ib1, ib2
c$$$          pos=ssite(ib)%pos
c$$$          ssite(ib)%pos0=pos
c$$$      enddo
c$$$      end subroutine cppos


c$$$      subroutine grdepl(nvrelx,indrlx,alpha,etot,irlxsh,pletot,dist)
c$$$C     implicit none
c$$$      integer irlxsh,nvrelx,indrlx(nvrelx)
c$$$      double precision pletot(6,2),etot,dist(9),alpha
c$$$      double precision grad,vec1(6)
c$$$      integer iv,ipm,ipv
c$$$      call info5(30,1,0,' GRDEPL: point %i of %i for grad shear: '//
c$$$     .     'etot=%d',irlxsh,2*nvrelx,etot,0,0)
c$$$C     Get index for current shear and store energy for that shear
c$$$      if (irlxsh .gt. 0) then
c$$$        iv = (irlxsh-1)/2 + 1
c$$$        ipv = iv
c$$$C     ipv = indrlx(iv)
c$$$        ipm = mod((irlxsh-1),2) + 1
c$$$        pletot(ipv,ipm) = etot
c$$$      endif
c$$$C     If this is last point, form gradient and exit
c$$$      if (irlxsh .eq. 2*nvrelx) then
c$$$        do  iv = 1, nvrelx
c$$$C     ipv = indrlx(iv)
c$$$          ipv = iv
c$$$          grad = (pletot(ipv,1) - pletot(ipv,2))/(2*alpha)
c$$$          pletot(ipv,1) = grad
c$$$        enddo
c$$$        irlxsh = 0
c$$$        return
c$$$      endif
c$$$C     Get shear index for next shear and whether + or -
c$$$      irlxsh = irlxsh+1
c$$$      iv = (irlxsh-1)/2 + 1
c$$$      ipv = indrlx(iv)
c$$$      ipm = mod((irlxsh-1),2) + 1
c$$$      if (ipv .lt. 1 .or. ipv .gt. 6)
c$$$     .     call rx('grdepl: something wrong with indrlx')
c$$$C     Make new shear
c$$$      call dvset(vec1,1,6,alpha)
c$$$      if (ipm .eq. 2) call dvset(vec1,1,6,-alpha)
c$$$      call dpzero(dist,6)
c$$$      call grdep2(iv,iv,indrlx,vec1,dist)
c$$$      dist(7) = 1
c$$$      end subroutine grdepl



