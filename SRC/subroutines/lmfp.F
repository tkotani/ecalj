      subroutine lmfp(llmfgw) 
      use m_lmfinit,only: irs1,irs2,lhf,bz_nevmx,bz_lmet,bz_w, irs1x10, 
     &     ctrl_lfrce,ctrl_mdprm,ctrl_maxit,ctrl_nbas,ctrl_nl,ctrl_nspec,ctrl_nspin,
     &     ctrl_zbak,ctrl_lrel,ctrl_ldos,ctrl_nitmv,ctrl_nvario,ctrl_quit,ctrl_tol,
     &     ham_ehf,ham_ehk,ham_seref,ham_eterms,
     &     lat_gam,lat_alat, lat_ldist,
     &     mix_umix,mix_tolu,mix_b,mix_bv,mix_w,mix_wc,mix_nsave,mix_mmix,
     &     sspec=>v_sspec, ssite=>v_ssite, 
     &     sstrnsymg, symgaf, nat


      use m_supot,only:sv_a_oorhat,zv_a_osmrho,Supot
      use m_mksym,only: rv_a_oag,rv_a_osymgr,iv_a_oistab,lat_nsgrp
      use m_suham,only: Suham
      use m_lattic,only: lat_dist,lat_plat,rv_a_opos
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
!!  aug2020. T.Kotani remove types s_* except s_spec and s_pot. Data are in m_lmfinit.
!       In future, we need clear protection after data are generated.
!! --> Old comments are removed to avoid confusions. See branches at github 
      implicit none
      integer ::iwdummy ,idummy
      integer procid,master,mpipid,nproc
      logical mlog
      logical lgors,cmdopt,bittst
      integer i,ifi,igets,iors,ipr,irs1x,isw,ix(5),j, !,iscr
     .     k,lcgf,leks,lgunit,lpnu,lrel,lrout,nbas,
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
C ... for relaxation
      logical ::xyzfrz(3),ltmp !,lshr=.false.
      integer:: icom , natrlx , nvrelx , ltb , itrlx , nm , irlxsh 
     .     , nitrlx , bitor , ng
      integer ,allocatable :: indrx_iv(:)
      real(8) ,allocatable :: w_rv(:)
c ki, for rlse and allocate oindrx
      integer,allocatable:: iv_tmp(:)
      double precision mdprm(6),gam(4),gam1,bstim,rhosig,pletot(6,2),
     .     plat0(3,3),dist0(9),dist(9)
      parameter (nm=3)
C  ... for LDA+U
      integer nlibu,lmaxu
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
!! Total energy mode --etot ==>moved to m_lmfinit ---
      call getpr(ipr)
      nbas = ctrl_nbas
      nl   = ctrl_nl
      nspec= ctrl_nspec
      nsp  = ctrl_nspin
      qbg = ctrl_zbak(1)
      maxit = int(ctrl_maxit)
      lrel =  ctrl_lrel/=0
      call setcc(lrel) !lrel=T means scaler relativistiv c=274.074d0 in a.u.
      if (lhf) maxit = 1
      nbaspp = nbas      
      stdo   = lgunit(1)
      irs1x = irs1
      plat = lat_plat
      if(ipr>=30) call praugm(sspec,0) !Printout properties of species
      
!!    --rs=3 => always read from atom file !2020aug.
!! Sep2020 takao: --rs=3 mode is removed. (--rs=3 meand fixed density Harris-foukner MD).
c      if ( iand(3,int(ctrl_lrs)) .eq. 3 ) irs1x = 0 
C     ... Setup for no screening transformation
!!  Sep2020 comment " Shorten site positions" here removed.
      
!! ... Setup for charge mixing
      dmxp= 0d0
      dmxp(2)=mix_b
      dmxp(9)=mix_bv
      dmxp(4:5)= mix_w(1:2)
      dmxp(3)= mix_wc
      dmxp(6 ) = int(mix_nsave)
      dmxp(7) = int(mix_mmix)
      call parms0(0,0,0d0,0)
C     ... Allocate memory for forces ! lfrce = int(ctrl_lfrce)
      if (ctrl_lfrce .ne. 0) then
        numq = 1
        if ( int(bz_lmet) .eq. 4 ) numq = 3
        allocate(frc_rv(3*nbas*numq))
      endif
C     ... Relaxation setup
      itrlx = 1
c      lat_ldist=0
      nitrlx=ctrl_nitmv
      mdprm=ctrl_mdprm
c      lshr = nint(mdprm(1)) .gt. 100 !lshr mode switch 
      if (nint(mdprm(1)) .eq. 0) nitrlx = 0
      if (nint(mdprm(1)) .gt. 0 .and. nint(mdprm(1)) .lt. 4) then
        call rx('lmf not set up for MD yet')
      endif
      if (nitrlx > 0) then !nitrlx >0 is for atomic position relaxiation
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
      endif
      
!!=== Re-entry for lshr (share mode); but share mode is currently commented out.
!!      Probably share mode is going to be moved to outside of lmfp =======
c    4 continue !
      if(procid==master) write(stdo,*)' lmfp : potential and basis setup ... '
      call supot() 
      call suham(sspec,ssite)
      if ( ctrl_quit == 8 ) call rx0('quit = ham')
      allocate(ftot_rv(3*nbas))
      nevmx = bz_nevmx !Maximum number of eigenvalues
      lrout = 1  ! Whether to evaluate output density and/or KS energy
      leks = 1
      j = 6
      if (bz_nevmx .eq. -1) then
        lrout = 0
        leks = 0
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
        call sudmtu ( nbas , nsp , nlibu , lmaxu , ssite , sspec , 0 !read dmatu from dmatu.ext
     .        , lldau_iv , ng , rv_a_osymgr , iv_a_oistab , dmatu , vorbdmat)
        if ( int(ctrl_quit) .eq. 16 ) call rx0('quit = dmat')
      endif
     
!! === Main iteration loops ===
      if(nitrlx>=1) allocate(pos_move(3,nbas))
      do 2000 itrlx = 1,max(1,nitrlx) ! for atomic position relaxiation (molecular dynamics) ===
         do 1000 iter = 1,max(1,maxit)! for electronic structure ===
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
               call rdovfa ( nbas, nspec, ssite, sspec,qbg, sv_a_oorhat )
!              read initial given by the superposition of spherical atomic densities by lmfa.
               nit1 = 0
            elseif ( 1<= mod(irs1x,4) .and. mod(irs1x,4)<=2 ) then 
               lbin = (irs1x==1) ! lbin=T:rst for --rs=1, or lbin=F:rsta for --rs=2,
               k = -1
               if (procid .eq. master) then
                  ifi=ifile_handle() 
                  if(lbin)       open(ifi,file='rst.'//trim(sname),form='unformatted')
                  if(.not. lbin) open(ifi,file='rsta.'//trim(sname)) 
               endif
               call mpibc1(ifi,1,2,mlog,'lmfp','ifi')
               bz_wx = bz_w
               k = iors ( ssite , sspec ,  fileid !read rst file
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
            if(ipr .ge. 50) then
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

C     Hang on to previous site density matrix for this iteration
            if (nlibu .gt. 0) then
               dmato=dmatu
               dmatu=0d0
            endif
            
c!! lda+u to recall previous +U potential stored in vorbdmat
c            inquire(file='vorbdmat.dat',exist=nexist,size=nsize)
c            if(iter==1.and.nexist.and.nsize>0) then
c               open(file='vorbdmat.dat',newunit=ifivorbdmat,form='unformatted')
c               read(ifivorbdmat) vorbdmat
c               close(ifivorbdmat)
c            endif

!! main part --- Make and diagonalize hamiltonian, make new charge density ---
            if(maxit==0) call info0(20,1,0,' lmfp  : zero iterations sought ... no band pass')
            call bndfp ( nlibu , lmaxu , lldau_iv , ssite , sspec ,    
     .           leks , lrout , ctrl_lfrce , lpnu , dmxp , iter , maxit  , !evl_rv, 
     .           frc_rv , dmatu, vorbdmat, llmfgw )
!! Check convergence of dmatu (density matrix for LDA+U) and update it. Get vorbdmat (U potential)
            if (nlibu .gt. 0 .and. maxit .gt. 0 .and. lrout .gt. 0) then
               call chkdmu ( nbas , nsp , nlibu , lmaxu , ssite , sspec !, sham 
     .              , 0 , dmatu , dmato, vorbdmat, tolu , umix
     .              , lldau_iv , ng , rv_a_osymgr , iv_a_oistab )
c               if(procid==master) then
c                 open(file='vorbdmat.dat',newunit=ifivorbdmat,form='unformatted')
c                 write(ifivorbdmat) vorbdmat
c                 close(ifivorbdmat)
c               endif
            endif
            
!! ==> --density (plot density) is moved to locpot.F(rho1mt and rho2mt) and mkpot.F(smooth part).
            
!! Write restart file (skip if --quit=band) ---
            if (procid .eq. master) then
               if ( ctrl_quit .ne. 4 ) then
!!  Suppress saving rst file in the middle of a shear (irlxsh > 0)
                  if (irs2.gt.0 .and. (lrout.gt.0 .or. maxit .eq. 0) .and. irlxsh .eq. 0) then
                     lbin = irs2 .ne. 2
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
            
c            if (cmdopt('--window=',9,0,strn)) call rx0('lmf : early exit (--window option)')
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
! kino's memo for call nwit              
!   lsc   :0 self-consistency achieved (diffe<=etol, qdiff=dmxp(11)<=qtol)
!         :1 if not self-consistent, but encountered max. no. iter.
!         :2 Harris energy from overlap of free atoms (iter=1 and lhf=t)
!         :3 otherwise
               call nwit ( int(ctrl_nvario), iter, maxit, lhf.or.irs1x==0.and.iter==1,
     &              leks+i , etol , qtol , dmxp(11) , 'cxhi' , amom , etot , lsc )
            endif
            call mpibc1(lsc,1,2,mlog,'lmfp','lsc')
            if (lsc .eq. 2 .and. .not. lhf .and. maxit .gt. 1) lsc = 3
            if (lsc .eq. 1 .and. lrout .gt. 0  .or. lsc .eq. 3) then
               if (iter .ge. maxit) lsc = 1
               if (iter .lt. maxit) lsc = 3
            endif
            if ( int(ctrl_quit) .eq. 4 ) call rx0 ( 'lmf : exit (--quit=band)' )
!!     Continue iterations toward self-consistency
            if (lsc .gt. 2) then
               irs1x=3
               cycle            !goto 5
            endif
            exit
 1000    continue ! ---------------- SCF (iteration) loop end ----
         if (nitrlx==0) exit
         
!! --- Molecular dynamics. not maintained recently but atomic position relaxation was working
!!         (last check is probably around 2015) ----------
         do ibas=1,nbas
            ssite(ibas)%pos0 = ssite(ibas)%pos
            pos_move(:,ibas) = ssite(ibas)%pos
         enddo
         mdprm = ctrl_mdprm
         
c$$$  else !without shear
!! main routine to relax atomic positions.  
         call relax ( ssite , sspec , itrlx , indrx_iv 
     .        , natrlx , nvrelx , frc_rv , p_rv , w_rv , 0 , 0d0 , pos_move, icom )
!!--- we removed a procedure at apr2021 2020
!!  Updated positions in site structure !WARN! t.kotani think this is confusing because
!!  CAUTION: 'positions written in ctrl' and 'positions written in rst' can be different.
         if(procid==master) write(ifipos) itrlx,pos_move
         do ibas=1,nbas
            ssite(ibas)%pos = pos_move(:,ibas)
         enddo
C     ... Exit when relaxation converged or maximum number of iterations
         if (icom .eq. 1) then !takao 
            if (procid .eq. master) then
               k = int(ctrl_nvario)
               flg = 'C67' !what? 
               call nwitsv(1+2,k,flg,nsp,amom,etot)
            endif
            if ( procid==master) then
               call tcx('lmfp')
               call fexit(0,111,' LMFP: relaxation converged after %i iteration(s)',itrlx)
            else
               call tcx('lmfp')
               call fexit(0,111,' ',0)
            endif
         endif
C     ... Restore minimum gradient positions if this is last step
         if (itrlx .eq. nitrlx) then
            call info0(20,1,0,' lmfp: restore positions for minimum g')
            call prelx1 ( 1 , nm , .false. , natrlx , nvrelx , indrx_iv , p_rv !lshr=F
     .           , pos_move )
C     Repack updated positions in site structure
            do ibas=1,nbas
               ssite(ibas)%pos = pos_move(:,ibas)
            enddo
         endif
!! New density after atom shifts 
!! If explicitly told to read from atom files after atom movmment
         irs1x=3
         call smshft(ssite,sspec,sv_a_oorhat,zv_a_osmrho,ctrl_lfrce) !sctrl,sham,slat,
!! Write restart file (to include new positions)
         if (procid .eq. master) then ! .and. .not. lshr) then
            ifi = ifile_handle()
            open(ifi,file='rst.'//trim(sname),form='unformatted') 
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
      if (allocated(lldau_iv)) deallocate(lldau_iv)
      if (allocated(ftot_rv)) deallocate(ftot_rv)
      if (allocated(p_rv)) deallocate(p_rv)
      if (allocated(w_rv)) deallocate(w_rv)
      if (allocated(indrx_iv)) deallocate(indrx_iv)
      if (allocated(frc_rv)) deallocate(frc_rv)
      call tcx('lmfp')
      end subroutine lmfp
