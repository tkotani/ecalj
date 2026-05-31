!> Calculate x0, \epsilon, spin susceptibility.
!!
!! eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
!! (only epsPP_lmf_chipm mode works).
module m_hx0fp0
  contains
subroutine hx0fp0()
  use m_ReadEfermi,only: Readefermi,ef
  use m_readqg,only:     Readqg,Readngmx2,ngpmx,ngcmx
  use m_hamindex,only:   Readhamindex
  use m_readeigen,only:  Readeval,Init_readeigen,Init_readeigen2
  use m_read_bzdata,only: Read_bzdata, nqbz,nqibz,qbz, wqt=>wt,q0i,nq0i,nq0ix,neps
  use m_genallcf_v3,only: Genallcf_v3,natom,nspin,nl,nn,nlnmx, nctot, alat, esmr, il,in,im,nlnm, plat, pos,ecore, tpioa
  use m_hamindex,only: ngrp
  use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
  use m_readqgcou,only: readqgcou
  use m_mpi,only: MPI__Initialize,MPI__root, &
       MPI__Broadcast,MPI__rank,MPI__size, MPI__consoleout,comm, &
     & MPI__InitQgroups, MPI__SplitXq, MPI__AutoSetup, &
     & comm_b => comm_b_xq, comm_k => comm_k_xq, &
     & mpi__root_k => mpi__root_k_xq, mpi__root_q, ipr, &
     & comm_q, comm_root_k => comm_root_k_xq, mpi__rank_root_k => mpi__rank_root_k_xq, &
     & iq_qgroup, n_qgroup, qgroup_root
  use m_rdpp,only: Rdpp, &   ! & NOTE: "call rdpp" generate following data.
       nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx,nprecx,mrecl,nblochpmx
  use m_zmel,only: Mptauof_zmel!, Setppovlz,Setppovlz_chipm   ! & NOTE: these data set are stored in this module, and used
  use m_itq,only: Setitq !set itq,ntq,nband,ngcmx,ngpmx to m_itq
  use m_freq,only: Getfreq3, getfreq2, &! & NOTE: call getfreq generate following data.
       frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw,niw !, frhis0,nwhis0 !output of getfreq
  use m_tetwt,only: Tetdeallocate,Gettetwt, &! & followings are output of 'L871:call gettetwt')
       whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb
  use m_w0w0i,only: W0w0i, w0,w0i ! w0 and w0i (head part at Gamma point)
  use m_wv_storage, only: wv_init_shm, wv_dump_shm_to_file, wv_dealloc
  use m_ll,only: ll
  use m_readgwinput,only: ReadGwinputKeys, ecut,ecuts,mtet,ebmx,nbmx,nmbas,imbas,egauss !nmbas is number of magnetic atoms
  use m_qbze,only: Setqbze, nqbze,nqibze,qbze,qibze
!  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_genallcf_v3,only: nprecb,mrecb,mrece,nqbzt,nband,mrecg
  use m_readVcoud,only: Readvcoud,vcousq,zcousq !,ngb,ngc
  use m_x0kf,only: x0kf_zxq,deallocatezxq,deallocatezxqi,zxqi,zxq
  use m_llw,only: WVRllwR,WVIllwI,MPI__sendllw2
  use m_w0w0i,only: w0w0i
  use m_lgunit,only:m_lgunit_init,stdo
  use m_readqg,only: Readqg0
!  use m_dpsion,only: dpsion5
  use m_gpu,only: gpu_init
  use m_ftox
  implicit none
  !! We calculate chi0 by the follwoing three steps.
  !!  gettetwt: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part
  !!  note: eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
  !!-------------------------------------------------
  ! cccc this may be wrong or correct cccccccccc
  !r Be careful for the indexing...
  !r      A routine idxlnmc(nindxv,nindxc,...  in index.f
  !r      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
  !r      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic).
  !r      The indexing starts with core first and then valence on top of core
  !r      So n-index in "in" for valence electron is different from "inv".
  real(8)    :: qp(3), quu(3), ua=1d0, vcmean, frr
  real(8)    :: schi=1d0, chg1, chg2, dumm1, dumm2
  logical    :: debug=.false., hx0, lqall
  logical    :: realomega=.true., imagomega=.true.
  logical    :: omitqbz=.false., chipm=.false., nolfco=.false., epsmode=.false., crpa=.false.
  integer    :: ixc, iqxini, iqxend, nwp, noccxv, ngb, ngc, ngrpx
  integer    :: i, is, iq, iw, ibas, imb, ibasx, lxx
  integer    :: ilmx, lb, nb, mb, ixx, ilm_r, nx_r
  integer    :: npr, ifif, ifwd, ifv, ierr
  integer    :: iqixc2, igb1, igb2, imb1, imb2
  integer    :: ifepsdatnolfc, ifepsdat, ifchipmn_mat
  integer    :: n_bpara, n_kpara, worker_inQtask, worker_auto
  character(11)       :: ttt
  character(128)      :: itag, outs=''
  character(3),  external :: charnum3
  logical,       external :: cmdopt0, cmdopt2
  integer,       external :: verbose
  real(8),    allocatable :: symope(:,:)
  integer,    allocatable :: nxx_r(:), aimbas(:)
  real(8),    allocatable :: svec(:,:), cvec(:,:), spinvec(:,:), consvec(:,:)
  real(8),    allocatable :: mmnorm(:), momsite(:)
  complex(8), allocatable :: zzr(:,:), epsi(:,:), x0mean(:,:,:)
  complex(8), allocatable :: epstinv(:,:), epstilde(:,:)
  call MPI__Initialize()
  call gpu_init(comm)
  call M_lgunit_init()
  call MPI__consoleout('hx0fp0')
  call cputid (0)
  if(verbose()>=100) debug= .TRUE.
  !! computational mode select ! takao keeps only the Sergey mode.
  if(ipr) write(stdo,"(a)") '--- Type numbers #1 #2 #3 [#2 and #3 are options] ---'
  if(ipr) write(stdo,"(a)") ' #1:run mode= 11: normal! 111: normal fullband! 10111 : normal  crpa!'
  if(ipr) write(stdo,"(a)") '             202: epsNoLFC! 203: eps!  222: chi^+- NoLFC'
  if(ipr) write(stdo,"(a)")  '-------------------------------------------------------'
  if(cmdopt2('--job=',outs)) then; read(outs,*) ixc
  elseif(MPI__root) then         ; read(5,*)    ixc
  endif
  call MPI__Broadcast(ixc)
  call cputid(0)
  !! List of Switches: !  normalm: normal eps mode; !  crpa: crpa mode !  epsmode: (normalm or epsmode)!  omitqbz: qbz>nqbz+1 are calulated
  !!  realomega: \chi on real axis  !  imagomega: \chi on imag omega !  lqall: limited range of \chi on real axis (mainly for memory reduction
  !!  chipm: \Chi_pm mode (nspin=2) !  nolfco: no local field correction
  lqall=.true.
  if(ixc==11) then;      if(ipr) write(stdo,*)"OK ixc=11  normal ";         epsmode=.false. ; lqall=.true.
  elseif(ixc==111) then; if(ipr) write(stdo,*)"OK ixc=111 normal fullband"; epsmode=.false.
  elseif(ixc==10011)then;if(ipr) write(stdo,*)"OK ixc=10011 crpa ";         epsmode=.false. ; crpa=.true.
  elseif(ixc==202) then; if(ipr) write(stdo,*)"OK ixc=202 eps NoLFC";       epsmode =.true. ; imagomega=.false.; omitqbz=.true.;nolfco=.true.
  elseif(ixc==203) then; if(ipr) write(stdo,*)"OK ixc=203 eps wLFC";        epsmode = .true.; imagomega=.false.; omitqbz=.true.
  elseif(ixc==222) then; if(ipr) write(stdo,*)"OK ixc=222 chipm noLFC";     epsmode = .true.; imagomega=.false.; omitqbz=.true.;nolfco=.true.
     chipm=.true.    !  elseif(ixc==12) realomega=.false.; ecorr_on=901; then ! Total energy test mode --> need fixing
  else; call rx( ' hx0fp0: given mode ixc is not appropriate')
  endif
  call Read_BZDATA(hx0)
  if(ipr) write(stdo,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp
  if(MPI__root.and.ipr) then
     do i=1,nqbz
        if(i<10 .OR. i>nqbz-10) write(stdo,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
        if(i==10 .AND. nqbz>18) write(stdo,"('... ')")
     enddo
     write(stdo,*)' nqbz nqibz =',nqbz,nqibz
  endif
  call Readefermi()
  if(ipr) write(stdo,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
  call genallcf_v3(incwfx=0) !use 'ForX0 for core' in GWIN
  if(chipm .AND. nspin==1) call rx( 'chipm mode is for nspin=2')
  if(nqbz /=nqbzt ) call rx(' hx0fp0_sc: nqbz /=nqbzt  in hbe.d')
  call ReadGWinputKeys()    !Readin GWinput
  !! Readin Offset Gamma --------  !      call ReadQ0P()
  !! Readin q+G. nqbze and nqibze are for adding Q0P related points to nqbz and nqibz.
  call Readngmx2() !return ngpmx and ngcmx in m_readqg
  call Setqbze()    ! extented BZ points list
  if(ipr) write(stdo,*)' num of zero weight q0p=',neps
  if(ipr) write(stdo,"(i3,f14.6,2x,3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
  if(ipr) write(stdo,"(' ngcmx ngpmx nqbz nq0i= ',2i8)") ngcmx,ngpmx,nqbz,nq0i
  !  do i = 1,nq0i+1; ini = nqbz*(i-1); do ix=1,nqbz;if(ipr) write(stdo,"('hx0fp0 qbze q0i=',i8,3f10.4,2x,3f10.4)") ini+ix,qbze(:,ini+ix);enddo
  !! Get space-group transformation information. See header of mptaouof.
  !! Here we use ngrpx=1 ==> "no symmetry operation in hx0fp0", c.f. hsfp0.sc.m.F case.
  !! ngrpx=1 (no symmetry operation in hx0fp0), whereas we use ngrp in eibzmode=T.
  ngrpx = 1
  allocate(symope(3,3),source=reshape([1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0],[3,3]))
  call Mptauof_zmel(symope,ngrpx) !no symmetry (ngrpx=1) for hx0fp0; use ngrp only in eibzmode
  call Setitq()
  iqxend = nqibz + nq0i
  write(stdo,ftox) ' nqibze nqibz nq0i=',nqibze,nqibz,nq0i
  call Readhamindex() ! Initialization of readEigen
  call init_readeigen() !EVU EVD are read in init_readeigen
  call init_readeigen2()
  if(verbose()>50) print *,'eeee exit of init_readeigen2'
  call Getfreq3(lqall,epsmode,realomega,imagomega,ua,mpi__root)
  writefreq_r: if(realomega .AND. mpi__root) then
     open(newunit=ifif,file='freq_r') !write number of frequency points nwp and frequensies in 'freq_r' file
     write(ifif,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
     do iw= nw_i,-1
        write(ifif,"(d23.15,2x,i6)") -freq_r(-iw),iw
     enddo
     do iw= 0,nw
        write(ifif,"(d23.15,2x,i6)") freq_r(iw),iw
     enddo
     close(ifif)
  endif writefreq_r
  if(ipr) write(stdo,"(' nw=',i5)") nw
  nwp = nw+1
  !! Get eigenvector corresponds to exp(iqr) (q is almost zero).
  if(epsmode) allocate(epsi(nw_i:nw,neps))
  Tetrahedroninitialization: block
    real(8):: ekt(nband,nqbze,nspin)
    do is = 1,nspin
       do iq = 1,nqbze
          ekt(:,iq,is)= readeval(qbze(:,iq),is)
       enddo
    enddo
    noccxv = maxval(count(ekt(1:nband,1:nqbze,1:nspin)<ef,1)) ! maximum no. occupied valence states
  endblock Tetrahedroninitialization
  if(noccxv>nband) call rx( 'hx0fp0: all the bands filled! too large Ef')
  if (MPI__root) then
     open(newunit=ifwd,file='__WV.d')
     write (ifwd,"(1x,10i14)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
     close(ifwd)
  endif
  if(omitqbz) then !! Set iqxini !omitqbz means skip loopf for iq=1,nqibz
     iqxini= nqibz + 1
  else
     iqxini= 1
  endif
  if(chipm ) then !transverse spin susceptibility
     allocate(aimbas(nmbas),source=abs(imbas(1:nmbas)))
     allocate( svec(nbloch,nmbas),source=0d0 )
     allocate( cvec(nbloch,nmbas),momsite(nmbas), mmnorm(nmbas)) !May2007
     cvec=0d0
     do imb=1,nmbas
        ibas= aimbas(imb)
        open (newunit=ifv,file='MixSpin.'//charnum3(ibas))
        read(ifv,*) ibasx,lxx
        allocate(nxx_r(0:lxx))
        do i=0,lxx
           read(ifv,*) nxx_r(i)   !   if(ipr) write(stdo,"(2i5,d13.6)") nxx_r(i)
        enddo
        allocate(spinvec((lxx+1)**2,maxval(nxx_r)),consvec((lxx+1)**2,maxval(nxx_r)))
        spinvec=0d0
        do ilmx = 1, (lxx+1)**2
           lb = ll(ilmx )         !  if(ipr) write(stdo,*)' lb=',lb,lxx,ilmx
           do ixx = 1, nxx_r(lb)  !  if(ipr) write(stdo,*)' nn=',nn,nxx_r(lb)
              if(ilmx==1) then
                 read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),chg1,chg2 ,consvec(ilmx,ixx)
              else
                 read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),dumm1,dumm2 ,consvec(ilmx,ixx)
              endif
           enddo
        enddo
        if(imb==1) schi=merge(1d0,-1d0, chg1-chg2>=0d0) !spin direction
        i=0
        if(ibas>1) i= sum(nblocha(1:ibas-1))
        do lb  = 0, lx (ibas) !!  ReOrdering of spinvec in natom ordering...
           do nb  = 1, nx (lb,ibas)
              do mb  = -lb, lb
                 i = i+1
                 ilmx = lb**2+ lb+ mb +1
                 svec(i,imb) = spinvec(ilmx,nb)
                 cvec(i,imb) = consvec(ilmx,nb)
                 if(ipr) write(stdo,"(' i lb mb svec svec**2=',3i4,2d13.5)") i,lb,mb,svec(i,imb),svec(i,imb)**2
              enddo
           enddo
        enddo
        deallocate(nxx_r,spinvec,consvec)
        close(ifv)
        mmnorm (imb) = sqrt(sum(svec(:,imb)**2))
        momsite(imb) = chg1-chg2
        if(ipr) write(stdo,"( 'mmom mmnorm= ',2f14.10)")  momsite(imb),mmnorm(imb)
     enddo
  endif

  ! ngb_max: max SHM size per q-point (nolfco→1, chipm→nmbas, normal→nblochpmx)
  ! nq_calc: number of q-points actually computed = iqxend-iqxini+1
  call MPI__AutoSetup(merge(1, merge(nmbas, nblochpmx, chipm), nolfco), &
                      nwhis, npm, niw, iqxend - iqxini + 1, &
                      worker_out=worker_auto, &
                      n_bpara_xq_out=n_bpara, n_kpara_xq_out=n_kpara)
  worker_inQtask = worker_auto
  if(ipr) write(stdo,'(1X,A,3I5)') 'MPI: worker_inQtask n_bpara n_kpara', worker_inQtask, n_bpara, n_kpara
  call MPI__InitQgroups(worker_inQtask)
  call MPI__SplitXq(n_bpara, n_kpara)  ! n_bpara>1: ω-parallel (comm_b splits freq range)
  if(ipr) write(stdo,ftox)'mpi_rank',mpi__rank,'iq_qgroup=',iq_qgroup,'n_qgroup=',n_qgroup
  if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0: sanity check. |q(iq=1)| /= 0')
  if(ipr) write(stdo,*)" chi_+- mode nolfc=",nolfco
  if(.NOT.chipm) allocate(zzr(1,1),source=(0d0,0d0)) !dummy
  iqloop: do iq = iqxini, iqxend  ! NOTE: qp=(0,0,0) is omitted when iqxini=2
    if( mod(iq - iqxini, n_qgroup) /= iq_qgroup ) cycle
    call cputid (0)
    qp  = qibze(:,iq)
    ! Readin diagonalized Coulomb interaction zcousq: E(\nu,I), Enu basis is given in PRB81,125102; vcousq: sqrt(v), as well.
    if(ipr) write(stdo,*); if(ipr) write(stdo,"('===== do 1001: iq qp=',i7,3f9.4,' ========')")iq,qp
    call Readqg0('QGcou',qp,   quu,ngc) ! ngc: the number of IPW for the interaction matrix (in QGcou),
    call Readvcoud(qp,iq,NoVcou=chipm) !Readin vcousq,zcousq ngb ngc for the Coulomb matrix
    ngb = ngc+nbloch
    if(ipr) write(stdo,"('  nbloch ngb ngc=',3i10)") nbloch,ngb,ngc
    if(chipm) then !npr is the dimension of zxq(npr,npr)
      npr = nmbas
    elseif(nolfco) then
      npr = 1
    else
      npr = ngb
    endif

    call wv_init_shm(npr, niw, (1-npm)*nwhis, nwhis, comm_q, mreclx=mrecl)
    if(epsmode) call writeepsopen()
    if(ipr) write(stdo,"(' ##### ',2i4,' out of nqibz+n0qi nsp=',2i4,' ##### ')")iq, nqibz + nq0i, nspin
    call x0kf_zxq(realomega,imagomega,qp,iq,npr,schi,crpa,chipm,nolfco,zzr,is_m_basis=.false.)
    if(mpi__root_k) then
      if(realomega) then
        if(     epsmode) call writerealeps()
        if(.NOT.epsmode) call WVRllwR(qp,iq,npr,npr,is_x0_m_basis=.false.,is_wc_m_basis=.false.)
        call deallocatezxq()
      endif
      if(imagomega) then
        if(     epsmode) call rx('hx0fp0: imagomega=T and epsmode=T is not implemented')
        if(.NOT.epsmode) call WVIllwI(qp,iq,npr,npr,is_x0_m_basis=.false.,is_wc_m_basis=.false.)
        call deallocatezxqi()
      endif
      if(.NOT.epsmode) then
        call MPI_barrier(comm_root_k, ierr)
        if (mpi__rank_root_k == 0) call wv_dump_shm_to_file(iq, mrecl, nblochpmx, realomega, imagomega)
      endif
    endif
    call mpi_barrier(comm_k, ierr)
    call mpi_barrier(comm_b, ierr)
  end do iqloop
  call wv_dealloc()
  call MPI_barrier(comm,ierr)
  if( .NOT. epsmode) call MPI__sendllw2(iqxend, n_qgroup, qgroup_root)
  !! == W(0) divergent part and W(0) non-analytic constant part.== Note that this is only for qp=0 -->iq=1
  !! get w0 and w0i (diagonal element at Gamma point.   !! This return w0, and w0i
  if(( .NOT. epsmode) .AND. MPI__rank==0) call w0w0i(nw_i,nw,nq0i,niw,q0i,is_wc_m_basis=.false.) !llw,llwI,
  ! === w0,w0i are stored to zw for qp=0 ===    !! === w_ks*wk are stored to zw for iq >nqibz ===
  call cputid(0)
  if(ixc==11)   call rx0( ' OK! hx0fp0 mode=11    read <Q0P> normal')
  if(ixc==111)  call rx0( ' OK! hx0fp0 mode=111   normal')
  if(ixc==10011)call rx0( ' OK! hx0fp0 mode=10011 crpa normal')
  if(ixc==202)  call rx0( ' OK! hx0fp0 mode=202   epsPP NoLFC')
  if(ixc==203)  call rx0( ' OK! hx0fp0 mode=203   eps LFC ')
  if(ixc==222)  call rx0( ' OK! hx0fp0 mode=222   chi+- NoLFC')
!  if(ixc==12)   call rx0( ' OK! hx0fp0 mode=12    Ecor mode')
contains
  subroutine writeepsopen()
    character(4), external :: charnum4
    itag=''
    if(cmdopt0('--interbandonly')) itag='.interbandonly'
    if(cmdopt0('--intrabandonly')) itag='.intrabandonly'
    iqixc2 = iq- (nqibz+nq0ix)
    if(( .NOT. chipm) .AND. nolfco) then
      if(allocated(x0mean)) deallocate(x0mean)
      allocate( x0mean(nw_i:nw,1,1) )
      x0mean=0d0
    endif
    if(mpi__root_q) then
      if(( .NOT. chipm) .AND. wqt(iq-nqibz)==0d0) then
        open(newunit=ifepsdatnolfc,file=trim('EPS'//charnum4(iqixc2)//'.nlfc.dat'//itag))
        write(ifepsdatnolfc,"(a)")'# qp(1:3)   w(Ry)   eps    epsi  --- NO LFC'
        if( .NOT. nolfco) then
          open(newunit=ifepsdat,file=trim('EPS'//charnum4(iqixc2)//'.dat'//itag))
          write(ifepsdat,"(a)") '# qp(1:3)   w(Ry)   eps  epsi --- LFC included. '
        endif
      endif
    endif
    if(chipm) then ! zzr is only for chipm.and.nolfco mode
      if( allocated(zzr)) deallocate(zzr,x0mean)
      allocate(zzr(ngb,nmbas),x0mean(nw_i:nw,nmbas,npr),source=(0d0,0d0))
      zzr(1:nbloch,1:nmbas) = svec(1:nbloch,1:nmbas)
    endif
    if(mpi__root_q) then
      if(chipm .AND. wqt(iq-nqibz)==0d0) then !! ... Open ChiPM* files for \Chi_+-
        open(newunit=ifchipmn_mat,file='ChiPM'//charnum4(iqixc2)//'.nlfc.mat')
        write(ifchipmn_mat,"(255i5)") nmbas
        write(ifchipmn_mat,"(255i5)") aimbas(1:nmbas)
        write(ifchipmn_mat,"(255e23.15)") momsite(1:nmbas)
        write(ifchipmn_mat,"(255e23.15)")  mmnorm(1:nmbas)
        write(ifchipmn_mat,"( ' Here was eiqrm: If needed, need to fix hx0fp0')")
      endif
    endif
  end subroutine writeepsopen
  subroutine writerealeps()
    use m_kind, only: kp => kindrcxq
    complex(kind=kp), allocatable :: x0meanx(:,:)
    !$acc update host (zxq)
    if(nolfco) forall(iw=nw_i:nw) x0mean(iw,:,:)=zxq(:,:,iw) !1x1
    if(nolfco .AND. ( .NOT. chipm)) then
      if (nspin==1) x0mean= 2d0*x0mean !if paramagnetic, multiply x0 by 2
      if (nspin==1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
    else
      if (nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
    endif
    if(nolfco) then
      ttt='without LFC'
    else
      ttt='with LFC'
    endif
    if(chipm) then
      if(ipr) write(stdo,*) '--- chi0_{+-}}^{-1}      --- '//ttt
    else
      if(ipr) write(stdo,*) '--- dielectric constant --- '//ttt
      if(ipr) write(stdo, *)" trace check for W-V"
    endif
    if(allocated(epstilde)) deallocate(epstilde,epstinv)
    allocate(epstilde(npr,npr),epstinv(npr,npr))
    iwloop: do iw = nw_i,nw
      frr= dsign(freq_r(abs(iw)),dble(iw))
      if( .NOT. chipm) then
        if(debug) write(stdo,*) 'xxx2 epsmode iq,iw=',iq,iw
        vcmean=vcousq(1)**2
        epsi(iw,iqixc2)= 1d0/(1d0 - vcmean*zxq(1,1,iw))
        if(mpi__root_q) then
          if(ipr) write(stdo,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2e23.15,3x, 2e23.15, &
               " vcmean x0mean =", 2e23.15,3x, 2e23.15)') iqixc2,iw,2*frr, &
               1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean, zxq(1,1,iw)
          write(ifepsdatnolfc,'(3f12.8,2x,e12.4,2e23.15,2x,2e23.15)') &
               qp, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
        endif
        if( .NOT. nolfco) then
          do igb1=1,npr
            do igb2=1,npr
              if(igb1==1 .AND. igb2==1) then
                epstilde(igb1,igb2)= -vcmean*zxq(igb1,igb2,iw)
              else
                epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
              endif
              if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
            enddo
          enddo
          epstinv(1:npr,1:npr)=epstilde(1:npr,1:npr)
          call matcinv(npr,epstinv(1:npr,1:npr))
          epsi(iw,iqixc2)= epstinv(1,1)
          if(mpi__root_q) then
            if(ipr) write(stdo,'( " iq iw omega eps epsi  wLFC=",2i6,f8.3,2e23.15,3x, 2e23.15)') &
                 iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
            if(ipr) write(stdo,*)
            write(ifepsdat,'(3f12.8,2x,e12.4,2e23.15,2x,2e23.15)') qp, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
          endif
        endif
      elseif(chipm) then ! ChiPM mode without LFC
        allocate( x0meanx(npr,npr) )
        x0meanx = cmplx(x0mean(iw,:,:), kind=kp) / 2d0 !in Ry unit.
        do imb1=1,npr
          do imb2=1,npr
            x0meanx(imb1,imb2) = x0meanx(imb1,imb2)/mmnorm(imb1)/mmnorm(imb2)
          enddo
        enddo
        if(mpi__root_q) write(ifchipmn_mat,'(3f12.8,2x,f20.15,2x,255e23.15)')qp, 2*schi*frr, x0meanx(:,:)
        deallocate(x0meanx)
      endif
    end do iwloop
    if(chipm) then
      close(ifchipmn_mat)
    else
      close(ifepsdatnolfc) ! = iclose( filepsnolfc)
      if( .NOT. nolfco) close(ifepsdat) !  = iclose(fileps)
    endif
  end subroutine writerealeps
endsubroutine hx0fp0
end module m_hx0fp0
