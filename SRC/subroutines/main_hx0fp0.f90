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
  use m_read_bzdata,only: Read_bzdata, ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,ginv, &
       dq_,qbz,wbz,qibz,wibz,qbzw, idtetf,ib1bz,idteti, nstar,irk,nstbz, wqt=>wt,q0i,nq0i ,nq0iadd,ixyz,nq0ix,neps
  use m_genallcf_v3,only: Genallcf_v3,natom,nspin,nl,nn,nlnmx, nctot, alat, esmr, il,in,im,nlnm, plat, pos,ecore, tpioa
  use m_hamindex,only: ngrp
  use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
  use m_readqgcou,only: readqgcou
  use m_mpi,only: MPI__Initialize,MPI__root, &
       MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size, MPI__consoleout,comm, &
     & MPI__SplitXq, MPI__Setnpr_col, comm_b, comm_k, mpi__root_k, mpi__root_q,ipr
  use m_rdpp,only: Rdpp, &   ! & NOTE: "call rdpp" generate following data.
       nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx,nprecx,mrecl,nblochpmx
  use m_zmel,only: Mptauof_zmel!, Setppovlz,Setppovlz_chipm   ! & NOTE: these data set are stored in this module, and used
  use m_itq,only: Setitq !set itq,ntq,nband,ngcmx,ngpmx to m_itq
  use m_freq,only: Getfreq3, getfreq2, &! & NOTE: call getfreq generate following data.
       frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw,niw !, frhis0,nwhis0 !output of getfreq
  use m_tetwt,only: Tetdeallocate,Gettetwt, &! & followings are output of 'L871:call gettetwt')
       whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb
  use m_w0w0i,only: W0w0i, w0,w0i ! w0 and w0i (head part at Gamma point)
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
  real(8):: qp(3),  qgbin(3),qx(3), ua=1d0 ! ua is a dummy.
  integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) ,icount,kold,isold, ndble=8, ngb,ngc !,nmbas
  real(8),allocatable:: vxcfp(:,:), wgt0(:,:)
  integer,allocatable :: ngvecpB(:,:,:),ngveccB(:,:), ngvecp(:,:), ngvecc(:,:)
  complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:), zw(:,:),zw0(:,:) 
  real(8),allocatable :: eqt(:), ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:), ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:)
  complex(8),allocatable :: trwv(:),trwv2(:) !,rcxq(:,:,:,:)
  complex(8) :: fff,img=(0d0,1d0)
  complex(8),allocatable :: wwk(:,:,:)
  integer,allocatable :: noccxvv(:) !n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:),
  real(8) ::qbzx(3),anfvec(3)
  logical :: debug=.false.
  integer,allocatable:: ibasf(:)
  real(8),allocatable :: transaf(:,:)
  logical :: realomega=.true., imagomega=.true.
  complex(8),allocatable:: epsi(:,:),gbvec(:),zzr(:,:),x0mean(:,:,:),zzr0(:)
  complex(8) :: epxxx,vcmean, vcmmmm
  complex(8),allocatable:: vcmmm(:)
  character(11) :: fileps,fileps23,filele
  character(16) :: filepsnolfc
  character(10) :: i2char
  character(20):: xxt
  real(8) :: Emin, Emax,emin2,emax2
  real(8) :: omg2max,omg1max,wemax
  real(8), allocatable :: freqr2(:)  , ekxxx(:,:,:)   !      logical::imagonly=.false.,realonly=.false. !,readgwinput
  integer::maxocc2, &
       ixc,iqxini,iqxend,iqxendx, &
                                !     &   ifhbe,    &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
       i,ngrpmx,mxx,ini,ix,ngrpx,&! & ngcmx,ngpmx nq0i,nq0ix,
       ndummy1,ndummy2,ifcphi,is,nwp, &! & ifvcfpout,,mdimx,nbloch
       ifepscond ,nw0,iw,ifinin,iw0,ifwwk,noccxv,noccx &
       ,ifwd,ifrcwi,ifrcw,ifianf,ibas &! & ,nprecx
       ,ibas1,irot,iq,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0dummy &
       ,kx,isf,kqxx,kp,job,noccxvx(2)=-9999,nwmax & ! & ,ifev1,ifev2 nbnbx,nhwtot,
       ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j   !     &   ,ngpmx !,  ifchipmlog
  real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm, dwry,dwh,omg_c,omg2
  integer:: incwfin,  verbose
  real(8):: quu(3), deltaq(3)!,qq(3) !,qqq(3)=0d0
  logical:: omitqbz=.false., noq0p
  logical,allocatable :: iwgt(:,:,:,:)
  complex(8),allocatable:: wgt(:,:,:)
  real(8),allocatable:: qbz2(:,:)
  logical :: qbzreg !if true, we use off-gamma mesh.
  integer,allocatable:: nstibz(:) !Nov2004 Miyake's tote
  real(8),allocatable:: ecqw(:,:) !,wiw(:)
  real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree,pi,efz,qfermi,alpha,rs,voltot,ecelgas,efx,valn
  integer:: iqbz,iqindx,iflegas,nmx,ifcor,nqitot,isx,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1
  real(8) :: wk4ec,faca !eclda_bh,eclda_pz,
  real(8),allocatable::    evall(:)
  complex(8),allocatable:: ovlpc(:,:),evecc(:,:)
  integer:: nev !,  ifdpin
  real(8),allocatable:: totexc(:), trpv(:),trlog(:) 
  integer:: necut,iecut
  integer:: ifv,lxx,ibasx,ilmx,ilm_r,nx_r,lb,nb,mb
  integer,allocatable:: nxx_r(:)
  real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
  character*3:: charnum3
  real(8)::chg1,chg2,spinmom,schi=1d0
  complex(8),allocatable:: ovlp(:,:),evec(:,:),ovlpi(:,:)
  real(8),allocatable::eval(:)
  integer:: new,nmxx,ii,iy,ipl1,ixx
  complex(8),allocatable :: ppovl(:,:),oo(:,:),x0inv(:,:),ppovlzinv(:,:)
  real(8)::qxx(3),ssm
  real(8),allocatable::SS(:),rwork(:),ss0(:)
  complex(8),allocatable:: UU(:,:),VT(:,:),work(:),zw0bk(:,:),ddd(:,:),vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:) 
  complex(8)::x0mx
  complex(8),allocatable:: UU0(:,:),VT0(:,:)
  logical ::  chipm=.false.,nolfco=.false.,epsmode=.false.,normalm=.false., crpa=.false. ,epsPPmode=.false.
  integer::  ife, idum4 
  real(8):: qs,qt,ww,muu, ddq(3)
  character(11) :: ttt
  integer:: nnmx,nomx
  ! Feb2006 time-reversal=off case
  logical :: timereversal, testtimer,onceww
  integer:: jpm,ncc
  real(8):: frr
  integer:: ipm,nrecoff
  real(8),allocatable:: ebb(:)
  logical :: evaltest 
  character*300:: aline
  integer:: imb,imb1,imb2!,nmbas_in 
  integer,allocatable:: aimbas(:), iibas(:) 
  complex(8),allocatable:: am1(:),am2(:),mmat(:,:), &
       x0mat(:,:),x0matinv(:,:),eiqrm(:)
  integer:: ifchipmn_mat !, ifchipm_fmat !,ifchipm_mat
  integer::ifstoner,ifx,i1
  real(8):: Istoner,zz1,zz2,zz3,zz4,Istoner0,jzero2,dumm1,dumm2
  complex(8):: trr,trr0,trr1     , zzzx(4,4), zzzy(4,4),trrx,mmatx(4,4),denom(4,4)
  real(8),allocatable:: eee(:),mmnorm(:), asvec(:,:),ssv(:,:),sproj(:,:),sprojx(:,:), momsite(:)
  real(8):: eex(4),eey(4),qvv(3)
  !      logical :: newaniso,newaniso2,newanisox !,z1offd
  integer :: ngb0,ifvcoud,idummy,ifepstinv,igb1,igb2,ngb_in,iq0,ifisk,iqx,ig,ifiss,iq0x
  complex(8),allocatable:: epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:)!zcousq(:,:),
  !      real(8),allocatable:: vcousq(:)
  real(8):: fourpi,sqfourpi,absq,vcou1,vcou1sq

  !! Eq.(40) in PRB81 125102
  !      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
  !     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
  complex(8),allocatable::sk(:),sks(:),skI(:),sksI(:), &
       w_k(:),w_ks(:),w_kI(:), w_ksI(:), s_vc(:),vw_k(:),vw_ks(:)
  integer:: lxklm,nlxklm,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,iwx,icc1,icc2!,niw,niwxxx,
  complex(8):: vc1vc2   !      integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
  integer,allocatable:: nwgt(:,:)
  real(8),allocatable:: aik(:,:,:,:)
  integer,allocatable:: aiktimer(:,:)
  integer:: l2nl
  logical:: tiii,iprintx=.false.,symmetrize,eibzmode !,eibz4x0
  real(8):: qread(3),imagweight,q00(3),rfac00,q1a,q2a
  character(128):: vcoudfile,aaax,itag
  integer:: src,dest
  logical:: lqall
  integer,allocatable ::  invgx(:) 
  integer:: k
  complex(8),allocatable:: ppovl_(:,:)
  logical:: readw0w0itest=.false.,hx0
  logical,external::cmdopt0
  integer:: ifq0p,ifwc,ifif,ierr,iqxx,ifi0,npr
  real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)
  logical:: cmdopt2
!,zmel0mode
  character(20):: outs=''
!  logical,save:: initzmel0=.true.
  real(8):: q0a,qa
!  complex(8),allocatable:: rcxq0(:,:,:,:)
  logical,allocatable::   mpi__task(:)
  integer,allocatable::   mpi__ranktab(:)
  integer :: n_kpara = 1, n_bpara = 1, npr_col, worker_inQtask, nqcalc
  call MPI__Initialize()
  call gpu_init(comm)
  call M_lgunit_init()
  call MPI__consoleout('hx0fp0')
  call cputid (0)
  if(verbose()>=100) debug= .TRUE. 
  hartree  = 2d0*rydberg()
  pi       = 4d0*datan(1d0)
  fourpi   = 4d0*pi
  sqfourpi = sqrt(fourpi)
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
  l2nl=2*(nl-1)
  allocate(symope(3,3),source=reshape([1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0],[3,3]))
  call Mptauof_zmel(symope,ngrpx) !we set thing in m_zmel for matrix elemenets generator
  !! ppbrd = radial integrals, cgr = rotated cg coeffecients (no rotatio here since nrgpx=1 for identity matrix)
  !! Rdpp gives ppbrd: radial integrals and cgr = rotated cg coeffecients.
  !!       --> call Rdpp(ngrpx,symope) is moved to Mptauof_zmel \in m_zmel
  !! Set nband, ngcmx, nqpmx and itq in m_zmel
  call Setitq()
  !! Pointer to optimal product basis
  !     nblochpmx = nbloch + ngcmx !rdpp \in mptauof_zmel \in m_zmel
  allocate(ngveccB(3,ngcmx)) 
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
  noccx  = noccxv + nctot
  if (MPI__root) then
     open(newunit=ifwd,file='__WV.d')
     write (ifwd,"(1x,10i14)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
     close(ifwd)
  endif
  allocate( zw(nblochpmx,nblochpmx) )
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

  n_bpara = 1
  if(cmdopt2('--nb=', outs)) read(outs,*) n_bpara
  if(nolfco) n_bpara = 1

  nqcalc = iqxend - iqxini + 1
!  if(cmdopt0('--zmel0')) nqcalc = nqcalc - 1
!  if(nqcalc < 1) call rx('hx0fp0: sanity check. nqcalc < 1: specify more than 2 q-points in zmel0 mode')
  n_kpara = max(mpi__size/(n_bpara*nqcalc), 1)  !Default setting of parallelization. b-parallel is 1.
  if(cmdopt2('--nk=', outs)) read(outs,*) n_kpara
  worker_inQtask = n_bpara * n_kpara
  if(ipr) write(stdo,'(1X,A,3I5)') 'MPI: worker_inQtask, n_bpara, n_kpara', worker_inQtask, n_bpara, n_kpara
  call MPI__SplitXq(n_bpara, n_kpara)

  allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz))
  allocate( nwgt(1,iqxini:iqxend))
  ! allocate( mpi__ranktab(iqxini:iqxend), source=[(mod(iq-1,mpi__size)           ,iq=iqxini,iqxend)])
  ! allocate( mpi__task(iqxini:iqxend),    source=[(mod(iq-1,mpi__size)==mpi__rank,iq=iqxini,iqxend)])
  allocate( mpi__ranktab(iqxini:iqxend), source=[(mod(iq-1,mpi__size/worker_inQtask)*worker_inQtask           ,iq=iqxini,iqxend)])
  allocate( mpi__task(iqxini:iqxend),    source=[(mod(iq-1,mpi__size/worker_inQtask)==mpi__rank/worker_inQtask,iq=iqxini,iqxend)])
  if(ipr) write(stdo,ftox)'mpi_rank',mpi__rank,'mpi__Qtask=',mpi__task
  if(ipr) write(stdo,ftox) 'mpi_qrank', mpi__ranktab
  !! llw, and llwI are for L(omega) for Q0P in PRB81,125102
!  allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )
  if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0.sc: sanity check. |q(iqx)| /= 0')
  if(ipr) write(stdo,*)" chi_+- mode nolfc=",nolfco
  if(.NOT.chipm) allocate(zzr(1,1),source=(0d0,0d0)) !dummy
  iqloop: do 1001 iq = iqxini,iqxend  ! NOTE: qp=(0,0,0) is omitted when iqxini=2
!    if(cmdopt0('--zmel0').and.iq==iqxini) cycle
    if( .NOT. MPI__task(iq) ) cycle
    call cputid (0)
    qp  = qibze(:,iq)
    q00= qibze(:,iqxini)
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
    call MPI__Setnpr_col(npr, npr_col) ! set the npr_col : split of npr(column) for MPI color_b
    if(epsmode) call writeepsopen()
    if(ipr) write(stdo,"(' ##### ',2i4,' out of nqibz+n0qi nsp=',2i4,' ##### ')")iq, nqibz + nq0i, nspin
    call x0kf_zxq(realomega,imagomega,qp,iq,npr,schi,crpa,chipm,nolfco, q00,zzr)
    if(mpi__root_k) then
    realomegamode: if(realomega) then !===RealOmega === W-V: WVR and WVI. Wing elemments: llw, llwi LLWR,LLWI
      if(mpi__root_k) then
        if(     epsmode) call writerealeps() !write eps file and close
        if(.NOT.epsmode) call WVRllwR(qp,iq,npr,npr_col)
        call deallocatezxq()
      endif
    endif realomegamode
    imagomegamode: if(imagomega) then ! ImagOmega start ============================
      if(mpi__root_k) then
        if(     epsmode) call rx('hx0fp0: imagoemga=T and epsmod=T is not implemented')
        if(.NOT.epsmode) call WVIllwI(qp,iq,npr,npr_col)
        call deallocatezxqi()
      endif
    endif imagomegamode
    endif
    call mpi_barrier(comm_k, ierr)
    call mpi_barrier(comm_b, ierr)
1001 enddo iqloop
  call MPI_barrier(comm,ierr)
  if( .NOT. epsmode) call MPI__sendllw2(iqxend,MPI__ranktab) !!! mpi send LLW to root.
  !! == W(0) divergent part and W(0) non-analytic constant part.== Note that this is only for qp=0 -->iq=1
  !! get w0 and w0i (diagonal element at Gamma point.   !! This return w0, and w0i
  if(( .NOT. epsmode) .AND. MPI__rank==0) call w0w0i(nw_i,nw,nq0i,niw,q0i) !llw,llwI,
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
    character*4:: charnum4
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
        write(ifepsdatnolfc,"(a)")' qp(1:3)   w(Ry)   eps    epsi  --- NO LFC'
        if( .NOT. nolfco) then
          open(newunit=ifepsdat,file=trim('EPS'//charnum4(iqixc2)//'.dat'//itag))
          write(ifepsdat,"(a)") ' qp(1:3)   w(Ry)   eps  epsi --- LFC included. '
        endif
      endif
    endif
    if(chipm) then ! zzr is only for chipm.and.nolfco mode
      if( allocated(zzr)) deallocate(zzr,x0mean)
      allocate(zzr(ngb,nmbas),x0mean(nw_i:nw,nmbas,npr_col),source=(0d0,0d0))
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
#ifdef __MP
    use m_mpi, only: MPI__GatherXqw => MPI__GatherXqw_kind4
#else
    use m_mpi, only: MPI__GatherXqw
#endif
    use m_kind,only: kp => kindrcxq
    complex(kind=kp), allocatable :: zxqw(:,:), x0meanx(:,:)
    allocate (zxqw(npr, npr))
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
    iq0 = iq - nqibz
    if(allocated(epstilde)) deallocate(epstilde,epstinv)
    allocate(epstilde(npr,npr),epstinv(npr,npr))
    iwloop: do 1015 iw  = nw_i,nw
      frr= dsign(freq_r(abs(iw)),dble(iw))
      call MPI__GatherXqw(zxq(:,:,iw), zxqw, npr, npr_col)
      if( .NOT. chipm) then
        if(debug) write(stdo,*) 'xxx2 epsmode iq,iw=',iq,iw
        vcmean=vcousq(1)**2 !fourpi/sum(qp**2*tpioa**2) !aug2012
        ! epsi(iw,iqixc2)= 1d0/(1d0 - vcmean*zxq(1,1,iw))
        epsi(iw,iqixc2)= 1d0/(1d0 - vcmean*zxqw(1,1))
        if(mpi__root_q) then
          if(ipr) write(stdo,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2e23.15,3x, 2e23.15, &
               " vcmean x0mean =", 2e23.15,3x, 2e23.15)') iqixc2,iw,2*frr, &
               1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean, zxqw(1,1) !x0mean(iw,1,1)
          write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') &
               qp, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
        endif
        if( .NOT. nolfco) then
          ix=0
          do igb1=ix+1,npr
            do igb2=ix+1,npr
              if(igb1==1 .AND. igb2==1) then
                epstilde(igb1,igb2)= -vcmean*zxqw(igb1,igb2) !aug2012
              else
                epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
              endif
              if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
            enddo
          enddo
          epstinv(ix+1:npr,ix+1:npr)=epstilde(ix+1:npr,ix+1:npr)
          call matcinv(npr-ix,epstinv(ix+1:npr,ix+1:npr))
          epsi(iw,iqixc2)= epstinv(1,1)
          if(mpi__root_q) then
            if(ipr) write(stdo,'( " iq iw omega eps epsi  wLFC=",2i6,f8.3,2e23.15,3x, 2e23.15)') &
                 iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
            if(ipr) write(stdo,*)
            write(ifepsdat,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') qp, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
          endif
        endif
      elseif(chipm) then ! ChiPM mode without LFC
        allocate( x0meanx(npr,npr) )
        call MPI__GatherXqw(cmplx(x0mean(iw,:,:),kind=kp), x0meanx, npr, npr_col)
        x0meanx = x0meanx/2d0 !in Ry unit.
        do imb1=1,npr
          do imb2=1,npr
            x0meanx(imb1,imb2) = x0meanx(imb1,imb2)/mmnorm(imb1)/mmnorm(imb2)
          enddo
        enddo
        if(mpi__root_q) write(ifchipmn_mat,'(3f12.8,2x,f20.15,2x,255e23.15)')qp, 2*schi*frr, x0meanx(:,:)
        deallocate(x0meanx)
      endif
1015 enddo iwloop
    if( allocated(gbvec) ) deallocate(gbvec)
    if(chipm) then
      close(ifchipmn_mat) 
    else
      close(ifepsdatnolfc) ! = iclose( filepsnolfc)
      if( .NOT. nolfco) close(ifepsdat) !  = iclose(fileps)
    endif
  end subroutine writerealeps
endsubroutine hx0fp0
end module m_hx0fp0

  !$$$!! --- legas mode is not working now. Need fixing... voltot ntot are not given.
  !$$$      if(epsmode.and.legas) then
  !$$$        call rx( ' LEGAS mode is not maintained well. Need some fixing.')
  !$$$        voltot=0d0
  !$$$        ntot=0d0
  !$$$        if(ipr) write(stdo,*)' Find LEGAS. legas =',legas
  !$$$        iflegas = 2101
  !$$$        open (iflegas,file='LEGAS')
  !$$$        read(iflegas,*)rs
  !$$$        close(iflegas)
  !$$$        alpha  = (9*pi/4d0)**(1d0/3d0)
  !$$$        qfermi = alpha/rs
  !$$$        efx  = qfermi**2
  !$$$        valn = efx**1.5d0*voltot/3d0/pi**2
  !$$$        write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
  !$$$        write (6,*)'     Exact Fermi momentum  qf  =', qfermi
  !$$$        write (6,*)'     Exact Fermi energy    Ef  =', efx
  !$$$        do iq = iqxini,iqxend ! q=(0,0,0) is omitted!
  !$$$          if(iq<=nqibz) cycle
  !$$$          if(ipr) write(stdo,*)' iq=',iq
  !$$$          iqixc2 = iq- (nqibz+nq0ix)
  !$$$          filele ='EPSEG'//charnum4(iqixc2)//'.dat'
  !$$$          ife = iopen ( filele,1,3,0)
  !$$$          write(ife,"(a)")
  !$$$     &          ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
  !$$$          q = qibze(:,iq)
  !$$$          qt= sqrt(sum(qibze(1:,iq)**2))*2d0*pi/alat
  !$$$          qs= qt/qfermi
  !$$$          if(ipr) write(stdo,"(' qs qfermi=',2d13.5)"    ) qs,qfermi
  !$$$          if(ipr) write(stdo,"(' q-q^2/2 q+q^2=',2d13.5)") qs-qs**2/2d0,qs+qs**2/2d0
  !$$$          do iw  = nw_i,nw
  !$$$            ww  = freq_r(iw)
  !$$$            muu = ww/qfermi**2
  !$$$            if(     qs<2d0 .and. muu < qs-qs**2/2d0) then
  !$$$              x0mx= -img*qfermi/(4*pi*qs)*2*muu
  !$$$            elseif( qs<2d0 .and. muu < qs+qs**2/2d0) then
  !$$$              x0mx= -img*qfermi/(4*pi*qs)*( 1d0-(muu/qs-.5d0*qs)**2 )
  !$$$            else
  !$$$              x0mx=0d0
  !$$$            endif
  !$$$            vcmmmm= 4*pi/qt**2
  !$$$            epsi(iw,iqixc2) = 1d0/(1- vcmmmm * x0mx)
  !$$$c            epsi(iw,iqixc2) = 1d0/(1- vcmmm(iq) * x0meanx)
  !$$$            write(ife,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)')
  !$$$     &        q, 2*ww,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
  !$$$          enddo
  !$$$        enddo
  !$$$        if(ipr) write(stdo,*)' ----------legas end--------'
  !$$$      endif

  !$$$!! Write TEECOR ecorr_on mode
  !$$$      if(imagomega.and.ecorr_on>0) then
  !$$$        hartree=2d0*rydberg()
  !$$$        ifcor   = iopen('TEECORR2',1,-1,0) ! output files
  !$$$        do iecut=1,necut
  !$$$          if(ipr) write(stdo,"( ' RPA Ec =' 3f23.15,'   ecut ecuts (Ry)=',2d12.4)")
  !$$$     &   totexc(iecut)*hartree,trpv(iecut)*hartree, trlog(iecut)*hartree
  !$$$     &    ,ecut(iecut),ecuts(iecut)
  !$$$          write(ifcor,*) '============================'
  !$$$          write(ifcor,*) 'Correlation energy Erpa (eV)'
  !$$$          write(ifcor,*) '============================'
  !$$$          write(ifcor,*)' ### '
  !$$$          write(ifcor,"(5e23.15)")
  !$$$     &     totexc(iecut)*hartree,trpv(iecut)*hartree,trlog(iecut)*hartree
  !$$$     &    ,ecut(iecut),ecuts(iecut)
  !$$$        enddo
  !$$$!! output ecqw !    write(ifcor,*)'### ecqw(q,w) ###'
  !$$$        write(ifcor,*)' nqibz =',nqibz
  !$$$        write(ifcor,*)' nq0i  =',nq0i
  !$$$        write(ifcor,*)' niw   =',niw
  !$$$        write(ifcor,*)' --- See details of Ec in ecor.chk ---'
  !$$$C... Write electron gas correlation energy
  !$$$c$$$        legas = .false.
  !$$$c$$$        INQUIRE (FILE = 'LEGAS', EXIST = legas)
  !$$$c$$$        if(legas) then !!! test for electron gas case.
  !$$$c$$$          call rx( ' LEGAS mode is not maintained well. Need some fixing.')
  !$$$c$$$          voltot=0d0
  !$$$c$$$          ntot=0d0
  !$$$c$$$          if(ipr) write(stdo,*)' find LEGAS. legas =',legas
  !$$$c$$$          iflegas = 2101
  !$$$c$$$          open (iflegas,file='LEGAS')
  !$$$c$$$          read(iflegas,*)rs
  !$$$c$$$          close(iflegas)
  !$$$c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
  !$$$c$$$          qfermi = alpha/rs
  !$$$c$$$          efx  = qfermi**2
  !$$$c$$$          valn = efx**1.5d0*voltot/3d0/pi**2
  !$$$c$$$          write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
  !$$$c$$$          write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
  !$$$c$$$          write (6,*)' egas  Exact Fermi energy    Ef  =', efx
  !$$$c$$$          if(tetra) call rx( 'legas You have to give ef of  tetrahedron')
  !$$$c$$$          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
  !$$$c$$$          qfermi= dsqrt(efz)
  !$$$c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
  !$$$c$$$          rs    = alpha/qfermi
  !$$$c$$$          write (ifcor,*)' --- electron gas ---'
  !$$$c$$$          write (ifcor,*)' density parameter rs= ', rs
  !$$$c$$$          write (ifcor,*)' kf= ',qfermi
  !$$$c$$$          write (ifcor,*)' ### Barth-Hedin formula'
  !$$$c$$$          ecelgas = eclda_bh(rs) * hartree * ntot
  !$$$c$$$          write (ifcor,*)ecelgas
  !$$$c$$$          write (ifcor,*)' ### Perdew-Zunger formula'
  !$$$c$$$          ecelgas = eclda_pz(rs) * hartree * ntot
  !$$$c$$$          write (ifcor,*)ecelgas
  !$$$c$$$          write (ifcor,*)' ### Gell-Mann and Brueckner formula'
  !$$$c$$$          ecelgas = (-0.0311d0 * dlog(rs) -0.048d0) * hartree * ntot
  !$$$c$$$          write (ifcor,*)ecelgas
  !$$$c$$$        endif
  !$$$      endif

! !--------------------------------------------------------------------
! real*8 function eclda_bh(rs)
!   real(8) :: rs,cp,rp,z
!   cp       = 0.0504d0*0.5d0 ! 0.5 changes unit from Ry to Hartree
!   rp       = 30.d0
!   z        = rs / rp
!   eclda_bh = -cp * ( (1.d0+z**3)*dlog(1.d0+1.d0/z) &
!        + 0.5d0*z - z**2 - 0.33333333d0 )
! END function eclda_bh
! !--------------------------------------------------------------------
! real*8 function eclda_pz(rs)
!   real(8) :: rs
!   if (rs >= 1.d0) then
!      eclda_pz = -0.1423d0 / (1.d0 + 1.0529d0*dsqrt(rs) + 0.334d0*rs)
!   else
!      eclda_pz = -0.0480d0 + 0.0311d0*dlog(rs) - 0.0116d0 * rs &
!           + 0.0020d0*rs*dlog(rs)
!   endif
! END function eclda_pz
! !--------------------------------------------------------------------
! subroutine wecqw(ifcor, &
!      nqibz,nqbz,nq0i,nqitot,niw, &
!      wibz,wqt,wx,freqx,ecqw)

!   implicit double precision (a-h,o-z)
!   dimension   wibz(nqibz),wqt(nq0i),wx(niw), &
!        freqx(niw),ecqw(nqitot,niw)
!   integer:: ifcor,nqibz,nqbz,nq0i,nqitot,ip,ix,niw
!   real(8):: rydberg
!   write(ifcor,*)'### ecqw(q,w) ###'
!   write(ifcor,*)'nqibz =',nqibz
!   write(ifcor,*)'nq0i  =',nq0i
!   write(ifcor,*)'niw   =',niw
!   do ip = 2,nqitot
!      if (ip <= nqibz) then
!         wk = wibz(ip)*0.5d0 ! 0.5 for the normalization of wibz
!      else
!         !        wk = wqt(ip-nqibz)*wibz(1)*0.5d0 ! 0.5 for the normalization of wibz
!         wk = wqt(ip-nqibz)* 1d0/dble(nqbz)
!      endif
!      write(ifcor,*)'### iq,wq = ',ip,wk
!      sume=0d0
!      do ix = 1,niw
!         write(ifcor,*)freqx(ix),ecqw(ip,ix),wx(ix)
!         sume=sume+  wx(ix)/(freqx(ix)*freqx(ix)) * ecqw(ip,ix)
!      enddo
!      write(ifcor,*) '  sum ecqw*wx=', wk*sume*2d0*rydberg()
!      ! end of ip-loop
!   enddo
!   return
! end subroutine wecqw
