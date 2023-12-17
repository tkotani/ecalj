!!  Calculate x0, \epsilon, spin susceptibility.
!!
!! eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
!! (only epsPP_lmf_chipm mode works).
program hx0fp0
  use m_ReadEfermi,only: Readefermi,ef
  use m_readqg,only:     Readqg,Readngmx2,ngpmx,ngcmx
  use m_hamindex,only:   Readhamindex
  use m_readeigen,only:  Readeval,Init_readeigen,Init_readeigen2
  use m_read_bzdata,only: Read_bzdata, ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,ginv, &
       dq_,qbz,wbz,qibz,wibz,qbzw, idtetf,ib1bz,idteti, nstar,irk,nstbz, wqt=>wt,q0i,nq0i ,nq0iadd,ixyz,nq0ix,neps
  use m_genallcf_v3,only: Genallcf_v3, nclass,natom,nspin,nl,nn, nlmto,nlnmx, nctot, alat, esmr,clabl,iclass, il,in,im,nlnm, &
       plat, pos,ecore, tpioa
  use m_hamindex,only: ngrp
  use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
  use m_readqgcou,only: readqgcou
  use m_mpi,only: MPI__hx0fp0_rankdivider2,MPI__task,MPI__Initialize,MPI__Finalize,MPI__root, &
       MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size, MPI__ranktab,MPI__consoleout,MPI__barrier
  use m_rdpp,only: Rdpp, &   ! & NOTE: "call rdpp" generate following data.
       nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx,nprecx,mrecl,nblochpmx
  use m_zmel,only: Mptauof_zmel, Setppovlz,Setppovlz_chipm   ! & NOTE: these data set are stored in this module, and used
  use m_itq,only: Setitq !set itq,ntq,nband,ngcmx,ngpmx to m_itq
  use m_freq,only: Getfreq3, &! & NOTE: call getfreq generate following data.
       frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw,niw !, frhis0,nwhis0 !output of getfreq
  use m_tetwt,only: Tetdeallocate,Gettetwt, &! & followings are output of 'L871:call gettetwt')
       whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb
  use m_w0w0i,only: W0w0i, w0,w0i ! w0 and w0i (head part at Gamma point)
  use m_ll,only: ll
  use m_readgwinput,only: ReadGwinputKeys, ecut,ecuts,mtet,ebmx,nbmx,nmbas,imbas,egauss
  use m_qbze,only: Setqbze, nqbze,nqibze,qbze,qibze
  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_readVcoud,only: Readvcoud,vcousq,zcousq,ngb,ngc
  use m_x0kf,only: X0kf_v4hz, X0kf_v4hz_init,ncount,kc !X0kf_v4hz_symmetrize,
!  use m_eibz,only:Seteibz, nwgt,neibz,igx,igxt,eibzsym
  use m_llw,only: WVRllwR,WVIllwI,MPI__sendllw2
  use m_w0w0i,only: w0w0i
  use m_lgunit,only:m_lgunit_init
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
  ! cccccccccccccccccccccccccccccccccccccccccccccccc
  real(8):: q(3),  qgbin(3),qx(3), ua=1d0 ! ua is a dummy.
  integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) ,icount,kold,isold, ndble=8
  real(8),allocatable:: vxcfp(:,:), wgt0(:,:)
  integer,allocatable :: ngvecpB(:,:,:),ngveccB(:,:), ngvecp(:,:), ngvecc(:,:)
  complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:), zw(:,:),zw0(:,:), zxq(:,:,:),zxqi(:,:,:)
  real(8),allocatable :: eqt(:), ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:), ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:)
  complex(8),allocatable :: trwv(:),trwv2(:),rcxq(:,:,:,:)
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
  character(11) :: fileps
  character(11) :: fileps23
  character(16) :: filepsnolfc
  character(11) ::  filele
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
       ,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas &! & ,nprecx
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
  real(8),allocatable:: totexc(:), trpv(:),trlog(:) !ecut(:),ecuts(:) ,
  integer:: necut,iecut
  integer:: ifv,lxx,ibasx,ilmx,ilm_r,nx_r,lb,nb,mb
  integer,allocatable:: nxx_r(:)
  real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
  character*3:: charnum3
  character*4:: charnum4
  complex(8),allocatable:: jcoup(:,:), mcm(:,:,:)
  real(8)::chg1,chg2,spinmom,schi=1d0
  complex(8),allocatable:: ovlp(:,:),evec(:,:),ovlpi(:,:)
  real(8),allocatable::eval(:)
  integer:: new,nmxx,ii,iy,ipl1,ixx
  complex(8),allocatable :: ppovl(:,:),oo(:,:),x0meanx(:,:),x0inv(:,:),ppovlzinv(:,:)
  real(8)::qxx(3),ssm
  real(8),allocatable::SS(:),rwork(:),ss0(:)
  complex(8),allocatable:: UU(:,:),VT(:,:),work(:),zw0bk(:,:),ddd(:,:),vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:) !,sqovlp(:,:),sqovlpi(:,:)
  integer::lwork,info,imin,ifzxq
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
  integer:: imb,imb1,imb2,nmbas_in 
  integer,allocatable:: aimbas(:), iibas(:) 
  complex(8),allocatable:: am1(:),am2(:),mmat(:,:), &
       x0mat(:,:),x0matinv(:,:),eiqrm(:)
  integer:: ifchipmn_mat, ifchipm_fmat !,ifchipm_mat
  integer::ifstoner,ifx,i1
  real(8):: Istoner,zz1,zz2,zz3,zz4,Istoner0,jzero2,dumm1,dumm2
  complex(8):: trr,trr0,trr1     , zzzx(4,4), zzzy(4,4),trrx,mmatx(4,4),denom(4,4)
  real(8),allocatable:: eee(:),mmnorm(:), asvec(:,:),ssv(:,:),sproj(:,:),sprojx(:,:), momsite(:)
  real(8):: eex(4),eey(4),qvv(3)
  !      logical :: newaniso,newaniso2,newanisox !,z1offd
  integer :: ngb0,ifvcoud,idummy,ifepstinv,igb1,igb2,ngb_in,nmbas1,nmbas2,iq0,ifisk,iqx,ig,nmbas1x,ifiss,iq0x
  complex(8),allocatable:: epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:)!zcousq(:,:),
  !      real(8),allocatable:: vcousq(:)
  real(8):: fourpi,sqfourpi,absq,vcou1,vcou1sq

  !! Eq.(40) in PRB81 125102
  !      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
  !     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
  complex(8),allocatable::sk(:),sks(:),skI(:),sksI(:), &
       w_k(:),w_ks(:),w_kI(:), w_ksI(:), s_vc(:),vw_k(:),vw_ks(:)
  complex(8),allocatable:: llw(:,:), llwI(:,:),aaamat(:,:)
  integer:: lxklm,nlxklm,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,iwx,icc1,icc2!,niw,niwxxx,
  complex(8):: vc1vc2   !      integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
  integer,allocatable:: nwgt(:,:)
  real(8),allocatable:: aik(:,:,:,:)
  integer,allocatable:: aiktimer(:,:)
  integer:: l2nl
  logical:: tiii,iprintx,symmetrize,eibzmode !,eibz4x0
  real(8):: qread(3),imagweight,q00(3),rfac00,q1a,q2a
  character(128):: vcoudfile,aaax,itag
  integer:: src,dest
  logical:: lqall
  integer,allocatable ::  invgx(:) !iclasst(:),
  integer:: ificlass,k
  complex(8),allocatable:: ppovl_(:,:)
  logical:: readw0w0itest=.false.,hx0,cmdopt0
  integer:: ifq0p,ifwc,ifif,ierr,iqxx,ifi0
  real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)
  logical:: cmdopt2,zmel0mode
  character(20):: outs=''
  logical,save:: initzmel0=.true.
  real(8):: q0a,qa
  complex(8),allocatable:: rcxq0(:,:,:,:)
  call MPI__Initialize()
  call M_lgunit_init()
  call MPI__consoleout('hx0fp0')
  call cputid (0)
  if(verbose()>=100) debug= .TRUE. 
  hartree  = 2d0*rydberg()
  pi       = 4d0*datan(1d0)
  fourpi   = 4d0*pi
  sqfourpi = sqrt(fourpi)
  !! computational mode select ! takao keeps only the Sergey mode.
  write(6,"(a)") '--- Type numbers #1 #2 #3 [#2 and #3 are options] ---'
  write(6,"(a)") ' #1:run mode= 11: normal! 111: normal fullband! 10111 : normal  crpa!'
  write(6,"(a)") '             202: epsNoLFC! 203: eps!  222: chi^+- NoLFC'
  write(6,"(a)")  '-------------------------------------------------------'
  if(cmdopt2('--job=',outs)) then; read(outs,*) ixc
  elseif(MPI__root) then         ; read(5,*)    ixc; endif
  call MPI__Broadcast(ixc)
  call cputid(0)
  !! List of Switches: !  normalm: normal eps mode; !  crpa: crpa mode !  epsmode: (normalm or epsmode)!  omitqbz: qbz>nqbz+1 are calulated
  !!  realomega: \chi on real axis  !  imagomega: \chi on imag omega !  lqall: limited range of \chi on real axis (mainly for memory reduction
  !!  chipm: \Chi_pm mode (nspin=2) !  nolfco: no local field correction
  lqall=.true.
  if(ixc==11) then;      write(6,*)"OK ixc=11  normal ";         epsmode=.false. ; lqall=.false.
  elseif(ixc==111) then; write(6,*)"OK ixc=111 normal fullband"; epsmode=.false. 
  elseif(ixc==10011)then;write(6,*)"OK ixc=10011 crpa ";         epsmode=.false. ; crpa=.true.
  elseif(ixc==202) then; write(6,*)"OK ixc=202 eps NoLFC";       epsmode =.true. ; imagomega=.false.; omitqbz=.true.;nolfco=.true.
  elseif(ixc==203) then; write(6,*)"OK ixc=203 eps wLFC";        epsmode = .true.; imagomega=.false.; omitqbz=.true.   
  elseif(ixc==222) then; write(6,*)"OK ixc=222 chipm noLFC";     epsmode = .true.; imagomega=.false.; omitqbz=.true.;nolfco=.true.
     chipm=.true.    !  elseif(ixc==12) realomega=.false.; ecorr_on=901; then ! Total energy test mode --> need fixing 
  else; call rx( ' hx0fp0: given mode ixc is not appropriate') ; endif
  call Read_BZDATA(hx0)
  write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp
  if(MPI__root) then
     do i=1,nqbz
        if(i<10 .OR. i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
        if(i==10 .AND. nqbz>18) write(6,"('... ')")
     enddo
     write(6,*)' nqbz nqibz =',nqbz,nqibz
  endif
  call Readefermi()
  write(6,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
  call genallcf_v3(incwfx=0) !use 'ForX0 for core' in GWIN
  if(chipm .AND. nspin==1) call rx( 'chipm mode is for nspin=2')
  if(nclass /= natom) call rx( ' nclass /= natom ') !! WE ASSUME iclass(iatom)= iatom
  call Readhbe()
  if(nqbz /=nqbzt ) call rx(' hx0fp0_sc: nqbz /=nqbzt  in hbe.d')
  if(nlmto/=nlmtot) call rx('hx0fp0: nlmto/=nlmtot in hbe.d')
  call ReadGWinputKeys()    !Readin GWinput
  !! Readin Offset Gamma --------  !      call ReadQ0P()
  !! Readin q+G. nqbze and nqibze are for adding Q0P related points to nqbz and nqibz.
  call Readngmx2() !return ngpmx and ngcmx in m_readqg
  call Setqbze()    ! extented BZ points list
  write(6,*)' num of zero weight q0p=',neps
  write(6,"(i3,f14.6,2x,3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
  write(6,"(' ngcmx ngpmx nqbz nq0i= ',2i8)") ngcmx,ngpmx,nqbz,nq0i
!  do i = 1,nq0i+1; ini = nqbz*(i-1); do ix=1,nqbz;write(6,"('hx0fp0 qbze q0i=',i8,3f10.4,2x,3f10.4)") ini+ix,qbze(:,ini+ix);enddo
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
  write(6,*) ' nqibz nqibze=',nqibz,nqibze
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
  if(MPI__root) write(6,"(' nw=',i5)") nw
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
     open(newunit=ifwd,file='WV.d')
     write (ifwd,"(1x,10i14)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
     close(ifwd)
  endif
  allocate( zw(nblochpmx,nblochpmx) )
  nspinmx = nspin
  if(omitqbz) then !! Set iqxini !omitqbz means skip loopf for iq=1,nqibz
     iqxini= nqibz + 1
  else
     iqxini= 1
  endif
!  if(cmdopt0('--rcxq0')) iqxend=iqxini
  if( chipm ) then !transverse spin susceptibility
     allocate(aimbas(nmbas))
     aimbas(1:nmbas)   = abs(imbas(1:nmbas))
     allocate(jcoup(nw_i:nw,neps) )
     allocate( svec(nbloch,nmbas) )  !sep2006
     svec=0d0
     allocate( cvec(nbloch,nmbas),momsite(nmbas), mmnorm(nmbas)) !May2007
     cvec=0d0
     do imb=1,nmbas
        ibas= aimbas(imb)
        open (newunit=ifv,file='MixSpin.'//charnum3(ibas))
        read(ifv,*) ibasx,lxx
        allocate(nxx_r(0:lxx))
        do i=0,lxx
           read(ifv,*) nxx_r(i)   !   write(6,"(2i5,d13.6)") nxx_r(i)
        enddo
        allocate(spinvec((lxx+1)**2,maxval(nxx_r)))
        allocate(consvec((lxx+1)**2,maxval(nxx_r)))
        spinvec=0d0
        do ilmx = 1, (lxx+1)**2
           lb = ll(ilmx )         !  write(6,*)' lb=',lb,lxx,ilmx
           do ixx = 1, nxx_r(lb)  !  write(6,*)' nn=',nn,nxx_r(lb)
              if(ilmx==1) then
                 read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),chg1,chg2 ,consvec(ilmx,ixx)
              else
                 read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),dumm1,dumm2 ,consvec(ilmx,ixx)
              endif
              !           write(6,"(2i5,d13.6)") ilmx, ixx, spinvec(ilmx,ixx)
           enddo
        enddo
        !! Calculate ChiPM. So sign of omega should be correct.
        if(imb==1) then !determine spin direction with respect to ibas=aimbas(imb=1)
           spinmom=(chg1-chg2)
           schi=1d0
           if(spinmom<0d0) then
              schi  = -1d0   ! This affects to dpsion. Obtained results
              ! should be the same in both mode.
           endif
        endif
        !!  ReOrdering of spinvec in natom ordering...
        i=0
        if(ibas>1) i= sum(nblocha(1:ibas-1))
        do lb  = 0, lx (ibas)
           do nb  = 1, nx (lb,ibas)
              do mb  = -lb, lb
                 i = i+1
                 ilmx = lb**2+ lb+ mb +1
                 svec(i,imb) = spinvec(ilmx,nb)
                 cvec(i,imb) = consvec(ilmx,nb)
                 write(6,"(' i lb mb svec svec**2=',3i4,2d13.5)") &
                      i,lb,mb,svec(i,imb),svec(i,imb)**2
              enddo
           enddo
        enddo
        deallocate(nxx_r,spinvec,consvec)
        close(ifv)
        mmnorm (imb) = sqrt(sum(svec(:,imb)**2))
        momsite(imb) = chg1-chg2
        write(6,"( 'mmom mmnorm= ',2f14.10)")  momsite(imb),mmnorm(imb)
     enddo
  endif
  allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz))

  !! === Use of symmetry. EIBZ procedure PRB81,125102 ===
  !! EIBZ mode memo for nolfco (right?)
  !! If eibzmode=T, it is efficient but can slightly break crystal symmetry.(how much?)
  !! This is because band connectivity is judged by just from band ordering in tetrahedron weitht tetwt5.
  !!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
  !! EIBZ mode
  eibzmode = .false. !eibz4x0()
  allocate( nwgt(1,iqxini:iqxend))
 
  !call seteibz(iqxini,iqxend,iprintx)
  !! Calculate x0(q,iw) and W == main loop 1001 for iq.
  !! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
  !! Thus it is necessary to do iq=1 in advance to performom iq >nqibz.
  !! (or need to modify do 1001 loop).
  !! iq>nqibz for ixc=11 is not time-consuming (right???)
  call MPI__hx0fp0_rankdivider2(iqxini,iqxend)

  !! llw, and llwI are for L(omega) for Q0P in PRB81,125102
  allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )
  !! ======== Loop over iq ================================
  if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0.sc: sanity check. |q(iqx)| /= 0')
  iqloop: do 1001 iq = iqxini,iqxend  ! NOTE: q=(0,0,0) is omitted when iqxini=2
     if(cmdopt0('--zmel0').and.iq==iqxini) cycle
     if( .NOT. MPI__task(iq) ) cycle
     call cputid (0)
     q  = qibze(:,iq)
     q00= qibze(:,iqxini)
     !! Readin diagonalized Coulomb interaction ===
     !! zcousq: E(\nu,I), Enu basis is given in PRB81,125102; vcousq: sqrt(v), as well.
     write(6,*)
     write(6,"('===== do 1001: iq q=',i7,3f9.4,' ========')")iq,q
     call Readvcoud(q,iq,chipm) ! Set ngb,ngc and vcousq,zcousq !for the Coulomb matrix
     if(( .NOT. chipm)) then
        if(allocated(zzr)) deallocate(zzr)
        allocate(zzr(1,1)) !dummy
        zzr=0d0
     endif
     write(6,"('  nbloch ngb ngc=',3i10)") nbloch,ngb,ngc
     !!
     if(chipm .AND. nolfco) then
        nmbas_in = nmbas
     elseif(nolfco) then
        nmbas_in = 1
     else
        nmbas_in = ngb
     endif
     nmbas1 = nmbas_in
     nmbas2 = nmbas1

     !     ! Open output files for epsilon mode =====
     itag=''
     if(cmdopt0('--interbandonly')) itag='.interbandonly'
     if(cmdopt0('--intrabandonly')) itag='.intrabandonly'

     if(epsmode) then
        iqixc2 = iq- (nqibz+nq0ix)
        if(( .NOT. chipm) .AND. nolfco) then
           allocate( x0mean(nw_i:nw,1,1) )
           x0mean=0d0
        endif
        if(( .NOT. chipm) .AND. wqt(iq-nqibz)==0d0) then
           open(newunit=ifepsdatnolfc,file=trim('EPS'//charnum4(iqixc2)//'.nlfc.dat'//itag))
           write(ifepsdatnolfc,"(a)")' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
           if( .NOT. nolfco) then
              open(newunit=ifepsdat,file=trim('EPS'//charnum4(iqixc2)//'.dat'//itag))
              write(ifepsdat,"(a)") ' q(1:3)   w(Ry)   eps  epsi --- LFC included. '
           endif
        endif
        !! zzr is only for chipm.and.nolfco mode
        if(chipm .AND. nolfco) then
           allocate(zzr(ngb,nmbas),x0mean(nw_i:nw,nmbas,nmbas))
           x0mean=0d0
           zzr   =0d0
           zzr(1:nbloch,1:nmbas) = svec(1:nbloch,1:nmbas)
        endif
        !! ... Open ChiPM* files for \Chi_+-
        if(chipm .AND. wqt(iq-nqibz)==0d0) then
           open(newunit=ifchipmn_mat,file='ChiPM'//charnum4(iqixc2)//'.nlfc.mat')
           write(ifchipmn_mat,"(255i5)") nmbas
           write(ifchipmn_mat,"(255i5)") aimbas(1:nmbas)
           write(ifchipmn_mat,"(255e23.15)") momsite(1:nmbas)
           write(ifchipmn_mat,"(255e23.15)")  mmnorm(1:nmbas)
           write(ifchipmn_mat,"( ' Here was eiqrm: If needed, need to fix hx0fp0')")
           if( .NOT. nolfco) then
              open(newunit=ifchipm_fmat,file='ChiPM'//charnum4(iqixc2)//'.fmat')
              write(ifchipm_fmat) nbloch, natom,nmbas, iqxini,iqxend, nw_i,nw
              write(ifchipm_fmat) aimbas(1:nmbas),momsite(1:nmbas),mmnorm(1:nmbas)
              write(ifchipm_fmat) nblocha(1:natom),svec(1:nbloch,1:nmbas)
              write(ifchipm_fmat) zzr0(1:nbloch) !zzr(1:nbloch,1)
           endif
        endif
     endif
     if(epsmode .AND. nolfco) then
        continue
     else
        write(6,*) "rcxq alloc ngb nwhis npm ---",ngb,nwhis,npm
        allocate( rcxq(ngb,ngb,nwhis,npm) )
     endif
     !! Set ppovlz: zmelt conversion. Non-orthogonality and Enu(coulomb-diagonal basis) related.
     !!     ppovlz is used in get_zmelt2 in m_zmel (called in x0kf_v4h).
     if(chipm .AND. nolfco) then
        call setppovlz_chipm(zzr,nmbas1)
     elseif(nolfco .AND. nmbas1==1) then !for <e^iqr|x0|e^iqr>
        call Setppovlz(q,matz=.true.)
     else                     !may2013  this removes O^-1 factor from zmelt
        call Setppovlz(q,matz=.true.) !.not.eibzmode)
     endif
     !! rcxq: imaginary part after x0kf_v4h and symmetrization.
     !! zxq ans zxqi are the main output after Hilbert transformation
     if(nolfco) then
        if(allocated(rcxq)) deallocate(rcxq)
        if(allocated(zxq) ) deallocate(zxq)
        if(allocated(zxqi) ) deallocate(zxqi)
        allocate( rcxq(nmbas1,nmbas2,nwhis,npm) )
        allocate( zxq (nmbas1,nmbas2,nw_i:nw), zxqi (nmbas1,nmbas2,niw))
     else
        allocate( zw0(ngb,ngb), zxq (ngb,ngb,nw_i:nw), zxqi(ngb,ngb,niw) )
     endif
     zxq=0d0;  zxqi=0d0;  rcxq = 0d0
     kold=-999
     isold=-999
     isloop: do 1003 is = 1,nspinmx
        write(6,"(' ##### ',2i4,' out of nqibz+n0qi nsp=',2i4,' ##### ')")iq, is, nqibz + nq0i,nspin
        if(debug) write(6,*)' niw nw=',niw,nw
        !  chi(charge) or chi_+-(spin when chipm=T)
        isf=is
        if(chipm) then
           write(6,*)" chi_+- mode nolfc=",nolfco
           if(is==1) isf=2
           if(is==2) isf=1
           rcxq=0d0
        endif
        do kx = 1, nqbz
           ekxx1(1:nband, kx)  = readeval(qbz(:,kx),   is )
           ekxx2(1:nband, kx)  = readeval(q+qbz(:,kx), isf)
        enddo
        call gettetwt(q,iq,is,isf,ekxx1,ekxx2,nband)
        !! == x0kf_v4hz is the main routine to accumalte imaginary part of x0 into rcxq ==
        epsppmode=  epsmode.and.nolfco
        write(6,*)'epsppmode=',epsppmode
        ierr = x0kf_v4hz_init(0, q, is, isf, iq, nmbas_in,crpa)
        ierr = x0kf_v4hz_init(1, q, is, isf, iq, nmbas_in,crpa)
        call x0kf_v4hz(q,is,isf,iq,nmbas_in,rcxq,epsppmode,iqxini,q00=q00) !,eibzmode
        call tetdeallocate() !--> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
        !  rcxq is the accumulating variable for spins
        !!    Symmetrize and convert to Enu basis by dconjg(tranpsoce(zcousq)*rcxq8zcousq if eibzmode
        ! if(is==nspinmx .OR. chipm) then !Apr2015. TK think " .OR. chipm" is required for chipm mode
        !    ! ecause rcxq is calculated for each is, symmetrized and its contribution
        !    ! s added to zxq in dpsion5.
        !    call x0kf_v4hz_symmetrize(q,iq,nolfco,zzr,nmbas_in,chipm,eibzmode,eibzsym(:,:,iq),rcxq)
        !    !  crystal symmetry of rcxq is recovered for EIBZ mode.
        ! endif
        if(debug) write(6,"(a)") ' --- goto dpsion5 --- '
        if(is==nspinmx .OR. chipm) then
           ! if(cmdopt0('--rcxq0')) then
           !    open(newunit=ifi0,file='rcxq0',form='unformatted')
           !    write(ifi0)rcxq
           !    close(ifi0)
           !    goto 1001
           ! elseif(cmdopt0('--zmel0')) then
           !    if(initzmel0) then
           !       open(newunit=ifi0,file='rcxq0',form='unformatted')
           !       allocate(rcxq0,mold=rcxq)
           !       read(ifi0)rcxq0
           !       close(ifi0)
           !       initzmel0=.false.
           !    endif
           !    q0a=sum(q00**2)**.5
           !    qa=sum(q**2)**.5
           !    if(abs(q0a-qa)>1d-12) rcxq = qa**2/(qa**2-q0a**2)*(rcxq - rcxq0)
           ! endif
           write(6,"('  nmbas1,nmbas2=',2i10)") nmbas1,nmbas2
           call dpsion5(realomega, imagomega, &
                rcxq, nmbas1,nmbas2, zxq, zxqi, &
                chipm, schi,is,  ecut,ecuts)
           if(nolfco .AND. epsmode) then
              do iw=nw_i,nw
                 x0mean(iw,:,:)=zxq(:,:,iw)
              enddo
           endif
           write(6,*)' --- end of dpsion5 ----',sum(abs(zxq)),sum(abs(zxqi))
        endif
        continue  
1003 enddo isloop
     if(allocated(rcxq)) deallocate(rcxq)
     realomegamode: if(realomega .AND. ( .NOT. epsmode)) then ! ===  RealOmega === W-V: WVR and WVI. Wing elemments: llw, llwi LLWR, LLWI
        call WVRllwR(q,iq,zxq,nmbas1,nmbas2)
        deallocate(zxq)
     elseif(realomega .AND. epsmode) then
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
           write(6,*) '--- chi0_{+-}}^{-1}      --- '//ttt
        else
           write(6,*) '--- dielectric constant --- '//ttt
           write(6, *)" trace check for W-V"
        endif
        !     ! prepare for iq0.
        iq0 = iq - nqibz
        if(allocated(epstilde)) deallocate(epstilde,epstinv)
        allocate(epstilde(ngb,ngb),epstinv(ngb,ngb))
        !     ! === iw loop for real axiw ===
        iwloop: do 1015 iw  = nw_i,nw
           frr= dsign(freq_r(abs(iw)),dble(iw))
           if( .NOT. chipm) then
              if(debug)write(6,*) 'xxx2 epsmode iq,iw=',iq,iw
              !     ! there is difference of two vcmean below since we use (sligthy) screened Coulomb (screenfac() in switch.F)
              !     !   NOTE that we use vcoul with screening (screenfac() is used in hvccfp0.F
              vcmean=vcousq(1)**2 !fourpi/sum(q**2*tpioa**2) !aug2012
              epsi(iw,iqixc2)= 1d0/(1d0 - vcmean*zxq(1,1,iw))
              write(6,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2e23.15,3x, 2e23.15, &
                   " vcmean x0mean =", 2e23.15,3x, 2e23.15)') iqixc2,iw,2*frr, &
                   1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean, zxq(1,1,iw) !x0mean(iw,1,1)
              write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') &
                   q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
              if( .NOT. nolfco) then
                 ix=0
                 do igb1=ix+1,ngb
                    do igb2=ix+1,ngb
                       if(igb1==1 .AND. igb2==1) then
                          epstilde(igb1,igb2)= -vcmean*zxq(igb1,igb2,iw) !aug2012
                       else
                          epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                       endif
                       if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                    enddo
                 enddo
                 epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                 call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                 epsi(iw,iqixc2)= epstinv(1,1)
                 write(6,'( " iq iw omega eps epsi  wLFC=",2i6,f8.3,2e23.15,3x, 2e23.15)') &
                      iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                 write(6,*)
                 write(ifepsdat,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
              endif
           elseif(chipm) then ! ChiPM mode without LFC
              allocate( x0meanx(nmbas,nmbas) )
              if(nolfco) then 
                 !$$$  c! --- three lines below may work for test purpose for legas. But not sure.
                 !$$$  c       vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
                 !$$$  c       write(ifchipmn,'(3f12.8,2x,f8.5,2x,2e23.15)')
                 !$$$  c     & q, 2*schi*frr, 1d0-vcmean*2*x0mean(iw,1,1)  !4*pi*alat**2/sum(q**2)/4d0/pi**2*x0mean(iw)
                 x0meanx = x0mean(iw,:,:)/2d0 !in Ry unit.
              else ! ChiPM mode with LFC... NoLFC part
                 zxq(1:ngb,1:ngb,iw) = zxq(1:ngb,1:ngb,iw)/2d0 ! in Ry.
                 do imb1=1,nmbas
                    do imb2=1,nmbas
                       x0meanx(imb1,imb2)= &
                            sum( svec(1:nbloch,imb1)* &
                            matmul(zxq(1:nbloch,1:nbloch,iw),svec(1:nbloch,imb2)))
                    enddo
                 enddo                 !     x0meanx= <m|chi^+-(\omega)|m>/<m|m>**2
              endif
              do imb1=1,nmbas
                 do imb2=1,nmbas
                    x0meanx(imb1,imb2) = x0meanx(imb1,imb2)/mmnorm(imb1)/mmnorm(imb2)
                 enddo
              enddo
              write(ifchipmn_mat,'(3f12.8,2x,f20.15,2x,255e23.15)')q, 2*schi*frr, x0meanx(:,:)
              if( .NOT. nolfco) write(ifchipm_fmat) q, 2*schi*frr, zxq(1:nbloch,1:nbloch,iw)
              deallocate(x0meanx)
           endif
1015    enddo iwloop
        if( allocated(zzr)   ) deallocate(zzr)
        if( allocated(x0mean)) deallocate(x0mean)
        if( allocated(gbvec) ) deallocate(gbvec)
        if(chipm) then
           close(ifchipmn_mat) !=iclose('ChiPM'//charnum4(iqixc2)//'.nlfc.mat')
           if( .NOT. nolfco) then
              close(ifchipm_fmat) !=iclose( 'ChiPM'//charnum4(iqixc2)//'.fmat')
           endif
        else
           filepsnolfc ='EPS'//charnum4(iqixc2)//'.nolfc.dat'
           close(ifepsdatnolfc) ! = iclose( filepsnolfc)
           if( .NOT. nolfco) then
              fileps = 'EPS'//charnum4(iqixc2)//'.dat'
              close(ifepsdat) !  = iclose(fileps)
           endif
        endif
     endif realomegamode
     imagomegamode: if (imagomega .AND. ( .NOT. epsmode)) then ! ImagOmega start ============================
        call WVIllwI(q,iq,zxqi,nmbas1,nmbas2)
        deallocate(zxqi)
     elseif(imagomega .AND. epsmode) then
        call rx('hx0fp0: imagoemga=T and epsmod=T is not implemented')
     endif imagomegamode
     if(allocated(vcoul)) deallocate(vcoul)
     if(allocated(zw0)) deallocate(zw0)
     if(allocated(zxq )) deallocate(zxq)
     if(allocated(zxqi)) deallocate(zxqi)
     if ( .NOT. epsmode) then
        close(ifrcwi)
        close(ifrcw)
     endif
1001 enddo iqloop
  call MPI__barrier()
!  if(cmdopt0('--rcxq0')) call rx0('end of --rcxq0 mode to generete rcxq0')
  if( .NOT. epsmode) call MPI__sendllw2(iqxend) !!! mpi send LLW to root.
  !! == W(0) divergent part and W(0) non-analytic constant part.==
  !!   Note that this is only for q=0 -->iq=1
  !! get w0 and w0i (diagonal element at Gamma point
  !! This return w0, and w0i
  if(( .NOT. epsmode) .AND. MPI__rank==0) call w0w0i(nw_i,nw,nq0i,niw,q0i) !llw,llwI,
  ! === w0,w0i are stored to zw for q=0 ===    !! === w_ks*wk are stored to zw for iq >nqibz ===
  call cputid(0)
  if(ixc==11)   call rx0( ' OK! hx0fp0 mode=11     read <Q0P> normal sergeyv')
  if(ixc==111)  call rx0( ' OK! hx0fp0 mode=111    normal sergeyv')
  if(ixc==10011)call rx0( ' OK! hx0fp0 mode=10011  crpa normal sergeyv')
  if(ixc==12)   call rx0( ' OK! hx0fp0 mode=12  Ecor sergeyv mode')
  if(ixc==101)  call rx0( ' OK! hx0fp0 mode=101 Ecor ')
  if(ixc==202)  call rx0( ' OK! hx0fp0 mode=202 sergeyv epsPP NoLFC')
  if(ixc==203)  call rx0( ' OK! hx0fp0 mode=203 sergeyv eps LFC ')
  if(ixc==222)  call rx0( ' OK! hx0fp0 mode=222 chi+- NoLFC sergeyv')
END PROGRAM hx0fp0

  !$$$!! --- legas mode is not working now. Need fixing... voltot ntot are not given.
  !$$$      if(epsmode.and.legas) then
  !$$$        call rx( ' LEGAS mode is not maintained well. Need some fixing.')
  !$$$        voltot=0d0
  !$$$        ntot=0d0
  !$$$        write(6,*)' Find LEGAS. legas =',legas
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
  !$$$          write(6,*)' iq=',iq
  !$$$          iqixc2 = iq- (nqibz+nq0ix)
  !$$$          filele ='EPSEG'//charnum4(iqixc2)//'.dat'
  !$$$          ife = iopen ( filele,1,3,0)
  !$$$          write(ife,"(a)")
  !$$$     &          ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
  !$$$          q = qibze(:,iq)
  !$$$          qt= sqrt(sum(qibze(1:,iq)**2))*2d0*pi/alat
  !$$$          qs= qt/qfermi
  !$$$          write(6,"(' qs qfermi=',2d13.5)"    ) qs,qfermi
  !$$$          write(6,"(' q-q^2/2 q+q^2=',2d13.5)") qs-qs**2/2d0,qs+qs**2/2d0
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
  !$$$        write(6,*)' ----------legas end--------'
  !$$$      endif

  !$$$!! Write TEECOR ecorr_on mode
  !$$$      if(imagomega.and.ecorr_on>0) then
  !$$$        hartree=2d0*rydberg()
  !$$$        ifcor   = iopen('TEECORR2',1,-1,0) ! output files
  !$$$        do iecut=1,necut
  !$$$          write(6,"( ' RPA Ec =' 3f23.15,'   ecut ecuts (Ry)=',2d12.4)")
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
  !$$$c$$$          write(6,*)' find LEGAS. legas =',legas
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
!   call cputid(0)
!   !      call MPI__Finalize
!   if(ixc==11)   call rx0( ' OK! hx0fp0 mode=11     read <Q0P> normal sergeyv')
!   if(ixc==111)  call rx0( ' OK! hx0fp0 mode=111    normal sergeyv')
!   if(ixc==10011)call rx0( ' OK! hx0fp0 mode=10011  crpa normal sergeyv')
!   if(ixc==12)   call rx0( ' OK! hx0fp0 mode=12  Ecor sergeyv mode')
!   if(ixc==101)  call rx0( ' OK! hx0fp0 mode=101 Ecor ')
!   if(ixc==202)  call rx0( ' OK! hx0fp0 mode=202 sergeyv epsPP NoLFC')
!   if(ixc==203)  call rx0( ' OK! hx0fp0 mode=203 sergeyv eps LFC ')
!   if(ixc==222)  call rx0( ' OK! hx0fp0 mode=222 chi+- NoLFC sergeyv')
! END PROGRAM hx0fp0

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
! !---------------------------------------------------------------------
! subroutine getsqovlp(q,ngc,ngb,sqovlp)
!   !! == Get sqrt of ppovl ==
!   implicit none
!   real(8)::q(3)
!   integer:: ngc,ngb,nbloch,i,nmxx,ix,iy,nev
!   complex(8):: sqovlp(ngb,ngb)
!   complex(8),allocatable:: ooo(:,:),ppo(:,:),sqovlpi(:,:),ppovl(:,:)
!   complex(8),allocatable:: ovlp(:,:),evec(:,:)
!   real(8),allocatable:: eval(:)
!   nbloch = ngb-ngc
!   if(ngc==0) goto 888

!   allocate(ppovl(1:ngc,1:ngc))
!   call readppovl0(q,ngc,ppovl)
!   allocate(ooo(ngc,ngc),ppo(ngc,ngc),evec(ngc,ngc),eval(ngc))
!   ooo= 0d0
!   do ix=1,ngc
!      ooo(ix,ix)=1d0
!   enddo
!   ppo = ppovl
!   deallocate(ppovl)
!   nmxx = ngc
!   evec = 0d0
!   eval = 0d0
!   call diagcv(ooo, ppo, &
!        evec, ngc, eval, nmxx, 1d99, nev)
!   write(6,*)' diagcv overlap ngc nev=',ngc,nev
!   deallocate(ooo,ppo)

! 888 continue
!   sqovlp=0d0
!   do i=1,nbloch
!      sqovlp(i,i)=1d0
!   enddo
!   do i=1,ngc
!      if(eval(i)<0d0) then
!         call rx( 'getsqovlp:  eval(i) <0d0')
!      endif
!      do ix=1,ngc;  do iy=1,ngc
!         sqovlp(ix+nbloch,iy+nbloch)= &
!              sqovlp(ix+nbloch,iy+nbloch) &
!              + evec(ix,i)* sqrt(eval(i))* dconjg(evec(iy,i))
!      enddo;      enddo
!   enddo
!   if(allocated(evec)) deallocate(evec)
!   if(allocated(eval)) deallocate(eval)
!   write(6,*)' end of getsqovlp'
!   !         sqovlpi = sqovlp
!   !         call matcinv(ngb,sqovlp)     !  inverse
!   !         ovlpi=ovlp
!   !         deallocate(ppovl,ovlp)
! end subroutine getsqovlp

!--------------------------------------------------------------------
!      subroutine test_xxx(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
!      implicit none
!      integer:: nblochpmx,nbloch,ngb,iw,i,iq
!      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
!      real(8):: freqq
!      logical :: smbasis
!      character*(*)::tagname
!      trwv2 = 0d0
!      forall( i = 1:ngb)
!        trwv2 = trwv2 + zw(i,i)
!      end forall
!      end
!--------------------------------------------------------------------

! !--------------------------------------------------------------------
! subroutine diagno00(nbloch,wpvc,eval)
!   !! == ontain eigenvalue only for input complex matrix wpvc(nbloch,nbloch)
!   implicit none
!   integer:: nbloch,nmx,nev,i
!   complex(8),allocatable:: ovlpc(:,:),evecc(:,:),wpvcc(:,:)
!   real(8)::emx,eval(nbloch)
!   complex(8):: wpvc(nbloch,nbloch)
!   allocate( ovlpc(nbloch,nbloch),evecc(nbloch,nbloch),wpvcc(nbloch,nbloch))
!   wpvcc= wpvc
!   ovlpc= 0d0
!   do i=1,nbloch
!      ovlpc(i,i)=1d0
!   enddo
!   eval=0d0
!   nev  = nbloch
!   nmx  = nbloch
!   call diagcv(ovlpc,wpvcc, evecc, nbloch, eval, nmx, 1d99, nev)
!   deallocate(ovlpc,evecc,wpvcc)
! end subroutine diagno00

