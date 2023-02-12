module m_lmfinit ! All ititial data (except rst/atm data via iors/rdovfa) !TK extensively rewrote at 2023feb.
  !                to run lmf-MPIK,lmfa,lmchk are stored in this module
  !                We perform 'call m_lmfinit_init', which sets all initial data.
  !! m_lmfinit_init have three stages. Search the word 'Block' in the followings.
  !! At the bottom of this code, which may/maynot be a help to read this code
  use m_ext,only :sname        ! sname contains extension. foobar of ctrl.foobar
  use m_struc_def,only: s_spec ! spec structures.
  use m_MPItk,only: master_mpi
  use m_lgunit,only: stdo,stdl
  use m_density,only: pnuall,pnzall !These are set here! log-derivative of radial functions.
  implicit none
  type(s_spec),allocatable:: v_sspec(:) !nspec: number of species in the cell
  integer,parameter::  noutmx=48,NULLI=-99999,nkap0=3,mxspec=256,lstrn=10000
  integer,parameter::  n0=10,nppn=2, nrmx=1501,nlmx=64 ,n00=n0*nkap0
  real(8),parameter::  fpi=16d0*datan(1d0), y0=1d0/dsqrt(fpi), pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi)
  real(8),parameter::  NULLR =-99999, fs = 20.67098d0, degK = 6.3333d-6 ! defaults for MD
  logical,parameter::  T=.true., F=.false.
  integer,protected::  io_help=0, lat_nkqmx,nat, lxcf 
  character(lstrn),protected:: sstrnsymg
  character(256),protected:: header,symg=' ',   symgaf=' '!for Antiferro
  logical,protected :: ham_frzwf,ham_ewald
  integer,protected:: nspc,procid,nproc,master=0,nspx,& 
       maxit,gga,ftmesh(3),nmto=0,lrsigx=0,nsp=1,lrel=1,lso=0 
  real(8),protected:: pmin(n0)=0d0,pmax(n0)=0d0,&
       tolft,scaledsigma, ham_oveps,ham_scaledsigma
  real(8):: cc !speed of light
  integer,protected :: smalit, lstonr(3)=0,nl 
  logical,protected :: lhf,lcd4
  !! ... STRUC
  real(8),protected:: dlat,alat=NULLR,dalat=NULLR,vol,avw 
  integer,protected:: nbas=NULLI,nspec
  !! ... SPEC
  real(8),protected:: omax1(3)=0d0,omax2(3)=0d0,wsrmax=0d0,sclwsr=0d0,vmtz(mxspec)=-.5d0
  character*8,allocatable,protected:: slabl(:)
  integer,protected:: lmxbx=-1,lmxax,nkaph,nkapi
  logical,allocatable,protected:: cstrmx(:),frzwfa(:)
  integer,allocatable,protected:: lmxb(:),lmxa(:),idmod(:,:),idu(:,:),kmxt(:),lfoca(:),lmxl(:),nr(:),& 
       nmcore(:), nkapii(:),nkaphh(:)
  real(8),allocatable,protected:: rsmh1(:,:),rsmh2(:,:),eh1(:,:),eh2(:,:), &
       rs3(:),alpha(:,:), uh(:,:),jh(:,:), eh3(:),&
       qpol(:,:),stni(:), pnusp(:,:,:),qnu(:,:,:), pnuspdefault(:,:),qnudefault(:,:),qnudummy(:,:), &
       coreq(:,:), rg(:),rsma(:),rfoca(:), rmt(:),pzsp(:,:,:), amom(:,:),spec_a(:),z(:),eref(:)
  character*(8),allocatable,protected:: coreh(:)
  !! ... SITE
  integer,allocatable,protected :: ispec(:)
  character(8),protected:: alabl
  real(8),allocatable,protected :: pos(:,:) 
  integer,allocatable,protected :: ifrlx(:,:),ndelta(:) 
  real(8),allocatable,protected :: delta(:,:),mpole(:),dpole(:,:)
  integer,allocatable,protected ::iantiferro(:)
  !! ... BZ
  integer,protected:: bz_lshft(3)=0, bz_lmet,bz_n,bz_lmull,bz_fsmommethod
  real(8),protected:: bz_efmax,bz_zval,bz_fsmom,bz_semsh(10),zbak,bz_range=5d0,bz_dosmax
  logical,protected:: bz_tetrahedron 
  !! ... Ewald
  real(8),protected:: lat_as,lat_tol,lat_rpad=0d0
  integer,protected:: lat_nkdmx
  !! ... STR
  real(8),protected:: str_rmax=nullr
  integer,protected:: str_mxnbr
  !! ... Iteration, MIX
  character(128),protected :: iter_mix=' '
  real(8),protected:: etol,qtol 
  integer,protected:: iter_maxit=1
  integer,protected:: mix_nsave 
  real(8),protected:: mix_tolu,mix_umix
  !!
  integer,protected:: pwmode,ncutovl ,ndimx    !ncutovl is by takao. not in lm-7.0beta.npwpad,
  real(8),protected:: pwemax,oveps,pwemin=0d0!,delta_stabilize
  integer,allocatable,protected ::  iv_a_oips (:)   ,  lpz(:),lpzex(:),lhh(:,:)
  !! ClebshGordon coefficient (GW part use clebsh_t)
  real(8) , allocatable,protected :: rv_a_ocg (:), rv_a_ocy (:)
  integer, allocatable,protected  :: iv_a_oidxcg(:), iv_a_ojcg (:)
  !!
  logical,protected:: addinv
  integer,protected:: ham_pwmode,ham_nkaph,ham_nlibu, nlmax,mxorb,lfrce
  integer,protected:: bz_nevmx, ham_nbf,ham_lsig,bz_nabcin(3)=NULLI, bz_ndos,ldos
  real(8),protected:: ham_seref, bz_w,lat_platin(3,3),lat_alat,lat_avw,lat_tolft,lat_gmaxin
  real(8),protected:: lat_gam(1:4)=[0d0,0d0,1d0,1d0] 
!  integer,protected:: lekkl  
  integer,protected:: lmaxu,nlibu
  integer,allocatable,protected::lldau(:)
  logical,protected:: lpztail=.false.
  integer,protected:: leks,lrout,plbnd,  pot_nlma,  pot_nlml,ham_nspx 
  integer,protected:: nlmto !total number of MTOs 
  real(8),protected:: socaxis(3) !SOC
  integer,protected:: natrlx,pdim
  logical,protected:: xyzfrz(3)
  integer,allocatable:: indrx_iv(:,:) 
  !! sspec unprotected, but there are changed only by readin parts, iors/rdovfa (see lmfp.f90)
  integer,allocatable,public,target:: ltabx(:,:),ktabx(:,:),offlx(:,:),ndimxx(:),norbx(:)
  integer,allocatable,public,protected:: jma(:),jnlml(:)
  !! DYN! molecular dynamics section DYN (only relaxiation, 2022-6-22)
  !   lrlxr: 0 no relaxation or dynamics
  !          4 relax with conjugate gradients
  !          5 relax with variable metric
  !          6 relax with Broyden
  !   rdhessr: T read hessian matrix
  !   xtolr: relaxation x-tolerance
  !   gtolr: relaxation g-tolerance
  !   stepr: step length
  !   nkillr: Remove hessian after this many steps
  integer,protected:: lrlxr,nkillr ! 2023feb
  logical:: rdhessr
  real(8),protected:: xtolr,gtolr,stepr    ! 2023feb
  integer,protected:: nitrlx ! nitmv   max number of mol-dynamics iterations
  ! mixrho
  logical,protected:: readpnu,v0fix,pnufix
  integer,protected:: broyinit,nmixinit,killj
  real(8),protected:: wtinit(3),wc,betainit 
  logical,protected:: bexist

  integer,parameter:: recln=512
  integer,protected:: reclnr,nrecs,nrecs2
  character(:),allocatable:: recrd(:)
contains
  subroutine m_lmfinit_init(prgnam) ! All the initial data are set in module variables from ctrlp.*
!    use m_toksw,only:tksw
!    use m_gtv,only: gtv,gtv_setrcd,gtv_setio
    use m_gtv2,only: gtv2_setrcd,rval2
    use m_cmdpath,only:cmdpath
    use m_ftox
    ! Inputs.
    !   ctrl file ctrl.sname
    !   prgnam: name of main program
    ! Outputs
    !    All the module variables. Only v_sspec can be modifed from outside (at readin part of atomic results).
    ! We have three stages (stage 1, stage 2 , stage 3) in this routine. Search 'stage'.
    ! Following memo is not so completed yet.
    !   BZ*  : Brillouin Zone related
    !   HAM* :  Hamiltonian related
    !   v_sspec : SPEC data.
    !   SITE: site information
    !   slabl : vector of species labels (species<-class<-nbas)
    !   avw   : the average Wigner-Seitz radius
    !   lrel  :specifies type of Schrodinger equation
    !         :0 nonrelativistic Schrodinger equation
    !         :1 scalar relativistic Schrodinger equation
    !   lxcf  :specifies type of XC potential. 
    !         :1 for Ceperly-Alder
    !         :2 for Barth-Hedin (ASW fit)
    !         :103 for PBE
    !   mxorb :nkaph \times (maximum number of lm channels in any sphere).
    !         MTO is specified by (n,l,m). (n=1,2,3. n=1:EH1, n=2:EH2, n=3:PZ)
    !   nbas  :number of atoms in the basis
    !   nkaph :The maximum number of radial functions centered at
    !         :particular R and l channel used in the lmto basis. +1 when we have lo
    !   nl    :1+Maximum l-cutoff for augmentation
    !   nkaphh :The maximum number of "principal quantum" numbers.
    !   nsp   :1 if not spin-polarized; otherwise
    !   nspc  :1 default. :2 spin-off diagonal included.
    !   nspec :number of species
    !   stde  :standard error file
    !   stdl  :standard log file
    !   stdo  :standard output file
    ! ----------------------------------------------------------------------
    implicit none
    include "mpif.h" 
    integer,parameter:: maxp=3
    character,intent(in)::  prgnam*(*)
    character strn*(recln),strn2*(recln)
    integer:: i_spec
    character fileid*64
    logical :: lgors,cmdopt,ltmp,ioorbp,cmdopt0
    double precision :: dval,dglob,xx(n0*2),dgets !,ekap(6)
    integer :: i,is,iprint, &
         iprt,isw,ifi,ix(n0*nkap0),j,k,l,lfrzw,lrs,k1,k2,mpipid 
    character*(8),allocatable::clabl(:)
    integer,allocatable:: ipc(:),initc(:),ics(:)
    real(8),allocatable:: pnuspc(:,:,:),qnuc(:,:,:,:),pp(:,:,:,:),ves(:),zc(:)
    integer:: dvec1(3)=1, dvec2(3)=0
    double precision :: orbp(n0,2,nkap0)
!    integer :: ohave,oics,opnusp,opp,oqnu,osgw,osoptc,oves,owk !osordn,
    real(8):: pnuspx(20) ,temp33(9)
    integer:: nnn
    real(8):: seref
    integer:: ib 
    integer,allocatable:: wowk(:)
    logical:: isanrg,l_dummy_isanrg
    integer:: lmxcg,lmxcy,lnjcg,lnxcg,nlm
    integer::nout,nn,i0,ivec(10),iosite
    integer:: io_tim(2),verbos!io_iactive,
    character(256)::  a,outs
    logical::  mlog=.false.
    integer:: lmxbj,lmxaj,nlbj
    double precision :: vsn,vers,xv(2*n0),xvv(3)
    character(256*16) :: bigstr=' '
    integer:: it
    logical::  debug=.false.
    integer:: lp1,lpzi
    real(8):: xxx
    real(8)::  avwsr
    integer:: ii,sw
    real(8):: dasum!,dglob
    character(128) :: nm
    real(8):: d2,plat(3,3),rydberg,rr
    real(8),allocatable ::rv(:)
    integer:: levelinit=0
    integer:: lx,lxx
    character*256:: sss,ch
    logical:: sexist
    integer:: ibas,ierr,lc, iqnu=0
    integer:: ifzbak,nn1,nn2,nnx,lmxxx,nlaj,isp
    integer,allocatable :: iv_a_oips_bakup (:)
    integer :: ctrl_nspec_bakup,inumaf,iin,iout,ik,iprior,ibp1,indx,iposn,m,nvi,nvl,nn1xx,nn2xx
    logical :: ipr10,fullmesh,lzz,fileexist
    logical:: logarr(100)
    integer,allocatable:: idxdn(:,:,:)
    if(master_mpi.and.cmdopt0('--help')) then
       call lmhelp(prgnam)
       call rx0('end of help mode')
    endif
    procid = mpipid(1)
    nproc  = mpipid(0)
    debug = cmdopt0('--debug')
    if(prgnam == 'LMF')    write(stdo,*) 'm_lmfinit:program LMF'
    if(prgnam == 'LMFGWD') write(stdo,*) 'm_lmfinit:program LMFGWD'
    if(prgnam == 'LMFA') write(stdo,*)   'm_lmfinit:program LMFA'
    ConvertCtrl2Ctrlp: block
      if(master_mpi) then
         inquire(file='ctrl.'//trim(sname),exist=fileexist)
         if( .NOT. fileexist) call rx("No ctrl file found!! ctrl."//trim(sname))
         GetCtrlp: block !Get ctrlp file
           character(512):: aaa,cmdl,argv
           aaa=''
           do i = 1, iargc()
              call getarg( i, argv )
              aaa=trim(aaa)//' '//trim(argv) !command inputs are connected
           enddo
           open(newunit=ifi,file='save.'//trim(sname),position='append')
           write(ifi,"(a)") 'Start '//trim(prgnam)//trim(aaa)
           close(ifi)
           cmdl=trim(cmdpath)//'ctrl2ctrlp.py '//trim(aaa)//'<ctrl.'//trim(sname)//' >ctrlp.'//trim(sname)
           write(stdo,*)'cmdl for python=',trim(cmdl)
           call system(cmdl) !Main part 
         endblock GetCtrlp
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr)
    end block ConvertCtrl2Ctrlp
    ReadCtrlp: block !Readin ctrlp given by ctrl2ctrlp.py above
      open(newunit=ifi,file='ctrlp.'//trim(sname))
      read(ifi,*) nrecs,reclnr,nrecs2
      allocate(character(reclnr):: recrd(nrecs2))
      do i = 1, nrecs2
         read(ifi,"(a)")recrd(i)
      enddo
      close(ifi)
    endblock ReadCtrlp
    Stage1GetCatok: block !Read Category-Token from recrd by rval2
      logical:: cmdopt0,cmdopt2,isanrg,parmxp
      integer:: setprint0,iprint,isw,ncp,nmix,broy,n,n1,n2,n3
      real(8):: avwsr,dasum,rydberg,wt(3),beta
      character(8):: fnam,xn
      call gtv2_setrcd(recrd)
      call rval2('STRUC_NSPEC',rr=rr, nreq=1);  nspec=nint(rr)
      call rval2('STRUC_NBAS', rr=rr, nreq=1);  nbas=nint(rr)
      call rval2('HAM_NSPIN',  rr=rr, defa=[real(8):: 1]);  nsp=nint(rr)
      allocate(pnuall(n0,nsp,nbas),pnzall(n0,nsp,nbas))
      allocate(pnusp(n0,nsp,nspec),qnu(n0,nsp,nspec), pzsp(n0,nsp,nspec),amom(n0,nspec),idmod(n0,nspec), &
           rsmh1(n0,nspec),eh1(n0,nspec),rsmh2(n0,nspec),eh2(n0,nspec), &
           stni(nspec), rg(nspec),rsma(nspec),rfoca(nspec),&!rcfa(2,nspec), & !,rsmfa(nspec)
           rmt(nspec),  spec_a(nspec),z(nspec),nr(nspec),eref(nspec), &
           coreh(nspec),coreq(2,nspec), idxdn(n0,nkap0,nspec), idu(4,nspec),uh(4,nspec),jh(4,nspec), &
           cstrmx(nspec),frzwfa(nspec), kmxt(nspec),lfoca(nspec),lmxl(nspec),lmxa(nspec),&
           lmxb(nspec),nmcore(nspec),rs3(nspec),eh3(nspec),&
           lpz(nspec),lpzex(nspec),nkapii(nspec),nkaphh(nspec),slabl(nspec),&
           pos(3,nbas),ispec(nbas),ifrlx(3,nbas),iantiferro(nbas))
      idu=0
      uh=0d0
      jh=0d0
      rs3=0.5d0
      eh3=0.5d0
      pnusp=0d0
      pzsp=0d0
      qnu=0d0
      lpz=0
      lpzex=0
      cstrmx=F !Exclude this species when auto-resizing sphere radii
      nkapii=1
      nkapi = 1
      lpzi = 0
      rsmh1 = 0d0
      rsmh2 = 0d0
      eh1  = 0d0
      eh2 = 0d0
      idmod = 0
      rfoca = 0d0
      rg = 0d0
      call rval2('IO_VERBOS' , rr=rr, defa=[real(8)::  30]); verbos=nint(rr)
      call rval2('IO_TIM'    , rr=rr, defa=[real(8)::  1 ]); io_tim=nint(rr)
      call rval2('STRUC_ALAT', rr=rr, nout=n);  alat=rr
      call rval2('STRUC_DALAT',rr=rr, nout=n);  dalat=rr
      call rval2('STRUC_PLAT', rv=rv, nreq=9);  plat=reshape(rv,shape(plat))
      call rval2('OPTIONS_HF' ,rr=rr, defa=[real(8):: 0]);  lhf= nint(rr)==1 ! for non-self-consistent Harris'
      call rval2('HAM_REL',    rr=rr, defa=[real(8):: 1]);  lrel=nint(rr)
      call rval2('HAM_SO',     rr=rr, defa=[real(8):: 0]);  lso=nint(rr)
      call rval2('HAM_SOCAXIS',rv=rv, defa=[0d0,0d0,1d0]);  socaxis=rv
      !  In general cases (except 001), we nees spin-symmetirc radial functions because we
      !  use  <up|Lz|up> for the place of <up|Lz|dn>, for example.
      !  If we like to take into account spin-dependent radial functions,
      !     We need to calculate ohsozz%soffd, ohsopm%sdiag. See 'call gaugm' at the end of augmat.
      call rval2('HAM_GMAX',   rr=rr, defa=[0d0]);          lat_gmaxin=rr
      call rval2('HAM_FTMESH', rv=rv, defa=[0d0,0d0,0d0]);  ftmesh=rv
      call rval2('HAM_TOL',   rr=rr,  defa=[1d-6]); tolft=rr
      call rval2('HAM_FRZWF', rr=rr,  defa=[real(8):: 0]); ham_frzwf= nint(rr)==1 
      call rval2('HAM_XCFUN', rr=rr,  defa=[real(8):: 2]); lxcf=nint(rr)
      call rval2('HAM_FORCES',rr=rr,  defa=[real(8):: 0]); lfrce=nint(rr)
      call rval2('HAM_RDSIG', rr=rr,  defa=[real(8):: 1]); lrsigx=nint(rr)
      call rval2('HAM_ScaledSigma', rr=rr, defa=[1d0]   ); scaledsigma=rr
      call rval2('HAM_EWALD', rr=rr,  defa=[real(8):: 0]); ham_ewald= nint(rr)==1
      call rval2('HAM_OVEPS', rr=rr,  defa=[1d-7]);   oveps=rr
      call rval2('HAM_PWMODE', rr=rr, defa=[real(8):: 0]);  pwmode=nint(rr)
      call rval2('HAM_PWEMAX', rr=rr, defa=[real(8):: 0]);  pwemax=rr
      call rval2('HAM_READP',rr=rr, defa=[real(8):: 0]); readpnu= nint(rr)==1
      call rval2('HAM_V0FIX',rr=rr, defa=[real(8):: 0]); v0fix =  nint(rr)==1
      call rval2('HAM_PNUFIX',rr=rr,defa=[real(8):: 0]); pnufix=  nint(rr)==1
      avw = avwsr(plat,alat,vol,nbas)
      specloop: do j=1,nspec !SPEC_ATOM j is spec index. In SPEC category, we do j=j+1 after we find ATOM=xx. See ctrl2ctrlp.py
         call rval2('SPEC_ATOM@'//xn(j), ch=ch); slabl(j)=trim(adjustl(ch))
         call rval2('SPEC_Z@'//xn(j),    rr=rr); z(j)=rr
         call rval2('SPEC_R@'  //xn(j), rr=rr, nout=n1); if(n1==1) rmt(j)=rr
         call rval2('SPEC_R/W@' //xn(j),rr=rr, nout=n2); if(n2==1.and.n1==0) rmt(j)=rr*avw
         call rval2('SPEC_R/A@' //xn(j),rr=rr, nout=n3); if(n3==1.and.n2==0.and.n1==0) rmt(j)=rr*alat
         call rval2('SPEC_A@'   //xn(j),rr=rr, defa =[0d0]);        spec_a(j)=rr
         call rval2('SPEC_NR@'  //xn(j),rr=rr, defa=[real(8)::0]);  nr(j)=nint(rr)
         call rval2('SPEC_RSMH@'//xn(j),rv=rv,nout=nnx); rsmh1(1:nnx,j)=rv       ! 1st MTO sets
         nn1 = findloc([(rsmh1(i,j)>0d0,i=1,nnx)],value=.true.,dim=1,back=.true.)! 1st MTO sets
         call rval2('SPEC_EH@'//xn(j),rv=rv,nreq=nn1); eh1(1:nn1,j)=rv           ! 1st MTO sets 
         call rval2('SPEC_RSMH2@'//xn(j),rv=rv,nout=nnx); rsmh2(1:nnx,j)=rv      ! 2nd MTO sets 
         nn2 = findloc([(rsmh2(i,j)>0d0,i=1,nnx)],value=.true.,dim=1,back=.true.)! 2nd MTO sets 
         call rval2('SPEC_EH2@'//xn(j), rv=rv, nreq=nn2); eh2(1:nn2,j)=rv        ! 2nd MTO sets 
         if(nn2>0) nkapi=2
         if(nn2>0) nkapii(j)=2
         call rval2('SPEC_LMX@'//xn(j), rr=rr, defa=[real(8):: 999]);     lmxb(j)=min(nint(rr), max(nn1,nn2)-1)
         call rval2('SPEC_LMXA@'//xn(j),rr=rr, defa=[real(8):: lmxb(j)]); lmxa(j)=nint(rr); if(rmt(j)==0d0) lmxa(j)=-1
         call rval2('SPEC_LMXL@'//xn(j),rr=rr, defa=[real(8):: lmxa(j)]); lmxl(j)=nint(rr) !'lmax for which to accumulate rho,V in sphere'
         call rval2('SPEC_P@'//xn(j),   rv=rv,nout=n); pnusp(1:n,1,j)=rv
         call rval2('SPEC_Q@'//xn(j),   rv=rv,nout=n); qnu(1:n,1,j)=rv
         if(nsp==2) then
            call rval2('SPEC_MMOM@'//xn(j), rv=rv,nout=n); qnu(1:n,2,j)=rv
         endif   
         call rval2('SPEC_NMCORE@'//xn(j),rr=rr, defa=[real(8):: 0]); nmcore(j)=nint(rr)
         call rval2('SPEC_PZ@'//xn(j),rv=rv,nout=n);  pzsp(1:n,1,j)=rv
         i0=1; if(z(j)<=8) i0=0
         call rval2('SPEC_LFOCA@'//xn(j),rr=rr, defa=[real(8):: i0]); lfoca(j)=nint(rr) !lfoca=0 or 1 
         call rval2('SPEC_KMXA@'//xn(j),rr=rr, defa=[real(8)::  3]); kmxt(j) =nint(rr)
         call rval2('SPEC_RSMA@'//xn(j),rr=rr, defa=[0.4d0*rmt(j)]); rsma(j) =rr
         call rval2('SPEC_IDMOD@'//xn(j),rv=rv,nout=n); idmod(1:n,j)=rv  !  'idmod=0 floats P to band CG, 1 freezes P, 2 freezes enu'
         call rval2('SPEC_FRZWF@'//xn(j),rr=rr,defa=[real(8):: 0]); frzwfa(j)= nint(rr)==1
         call rval2('SPEC_IDU@'//xn(j), rv=rv,nout=n); idu(1:n,j)=nint(rv)!U mode: 0 nothing, 1 AMF, 2 FLL, 3 mixed. +10:no LDA+U if sigm.* exist
         ! 2019 automatic turn off lda+u mode; Enforce uh=jh=0 when sigm exist !!
         call rval2('SPEC_UH@'//xn(j),  rv=rv,nout=n); uh(1:n,j)=rv ! Hubbard U for LDA+U
         call rval2('SPEC_JH@'//xn(j),  rv=rv,nout=n); jh(1:n,j)=rv ! Exchange parameter J for LDA+U
         call rval2('SPEC_C-HOLE@'//xn(j),ch=ch); coreh(j)=trim(adjustl(ch)) ! Channel for core hole
         call rval2('SPEC_C-HQ@'//xn(j),rv=rv,defa=[-1d0,0d0]) ; coreq(:,j)=rv
!              'Charge and Moment of core hole:'// &
!               Q(spin1) = full + C-HQ(1)/2 + C-HQ(2)/2 
!               Q(spin2) = full + C-HQ(1)/2 - C-HQ(2)/2
         call rval2('SPEC_ATOM_EREF'//xn(j), rr=rr,defa=[0d0]); eref(j)=rr
      enddo specloop
      lmxbx=maxval(lmxb)
      ibasloop: do j = 1, nbas
         call rval2('SITE_ATOM@'//xn(j),ch=ch); alabl=trim(adjustl(ch))
         ispec(j)= findloc([(trim(alabl)==trim(slabl(i)),i=1,nspec)],dim=1,value=.true.)
         if(ispec(j)==0) call rx('Category SITE referred to nonexistent species: '//trim(alabl))
         call rval2('SITE_POS@'//xn(j),rv=rv,  nout=n); pos(1:n,j)=rv !cartesian in alat
         if(n/=3) then
            call rval2('SITE_XPOS@'//xn(j),rv=rv, nout=n) !fractional (POSCAR direct) in alat
            if(n/=3) call rx('SITE_POS is not supplied')
            pos(:,j)= matmul(plat,rv)
         endif   
         call rval2('SITE_RELAX@'//xn(j),rv=rv,defa=[real(8):: 1,1,1]); ifrlx(:,j)=nint(rv)!relax site positions (lattice dynamics) 
         call rval2('SITE_AF@'//xn(j),   rr=rr,defa=[real(8):: 0]); iantiferro(j)=nint(rr) !AF=i and AF=-i should be AFerro-pair.
         !                                                                                  We look for space-group operation with spin-flip')
      enddo ibasloop
      call rval2('STR_RMAXS',rr=rr,nout=n); str_rmax=rr  ! Radial cutoff for strux, in a.u.'
      if(n ==0) then
         call rval2('STR_RMAX',rr=rr)! 'Radial cutoff for strux, in units of avw
         str_rmax=rr*avw  
      endif
      call rval2('STR_MXNBR',rr=rr, defa=[real(8):: 0d0]); str_mxnbr=rr !'Max number of nbrs (for dimensioning arrays)')
      call rval2('BZ_NKABC',rv=rv, nout=n); bz_nabcin(1:n)=nint(rv) !'No. qp along each of 3 lattice vectors.'//new_line('a')//'   '//&
      call fill3in(n,bz_nabcin) !filled to the end if n<3
      call rval2('BZ_BZJOB',rv=rv, nout=n); bz_lshft(1:n)=nint(rv) !  '0 centers BZ mesh at origin, 1 centers off origin'// &
      call fill3in(n,bz_lshft)
      call rval2('BZ_METAL', rr=rr, defa=[real(8):: 3]); bz_lmet=nint(rr) !'0 insulator only; 3 for metal (2 is for maintenance)')
      call rval2('BZ_TETRA', rr=rr, defa=[real(8):: 1]); bz_tetrahedron= nint(rr)==1 ! & tetrahedron switch
      if(cmdopt0('--tdos') .OR. cmdopt0('--pdos')) then
         write(stdo,*)' --tdos or --pdos enforces BZ_METAL=3 and BZ_TETRA=T'
         bz_lmet=3
         bz_tetrahedron=.true.
      endif
      call rval2('BZ_N', rr=rr, defa=[real(8):: 0]); bz_n=nint(rr) 
      !          N>0: Polynomial order for Methfessel-Paxton sampling', N=0: Conventional Gaussian sampling'// &
      !          N<0: Broadening by Fermi-Dirac distribution'.     Use in conjunction with BZ_W
      call rval2('BZ_W', rr=rr, defa=[5d-3]); bz_w=rr
      ! For BZ_N>=0, Line broadening for sampling integratio. For BZ_N<0, Temperature for Fermi distribution (Ry)
      call rval2('BZ_ZBAK',  rr=rr, defa=[0d0]); zbak=rr !'Homogeneous background charge'
      call rval2('BZ_SAVDOS',rr=rr, defa=[real(8):: 0]); ldos=nint(rr)! '0(F) or 1(T): Write dos.tot.* file (settings are NPTS and DOS)'
      call rval2('BZ_NPTS',  rr=rr, defa=[real(8):: 2001]); bz_ndos=nint(rr) !'No. DOS points (sampling integration)')
      call rval2('BZ_DOSMAX',rr=rr, defa=[40d0/rydberg()]); bz_dosmax=rr ! Maximum energy to which DOS accumulated, relative to Efermi
      call rval2('BZ_EFMAX', rr=rr, defa=[5d0]); bz_efmax=rr !Find evecs up to efmax'
      call rval2('BZ_NEVMX', rr=rr, defa=[real(8):: 0]); bz_nevmx=nint(rr)! Find at most nevmx eigenvec. NEVMX=0:internal default, NEVMX<0:no eigenvecs
      if( cmdopt0('--tdos') .OR. cmdopt0('--pdos') .OR. cmdopt0('--zmel0')) bz_nevmx=999999
      call rval2('BZ_FSMOM',rr=rr, defa=[NULLR]); bz_fsmom=rr !'Fixed-spin moment (fixed-spin moment method)')
      call rval2('BZ_FSMOMMETHOD', rr=rr, defa=[real(8):: 0]); bz_fsmommethod=nint(rr) !'Method of Fixed-spin moment 0:original 1:discrete')
      call rval2('SYMGRP', ch=ch); symg=adjustl(ch) ! Generators for symmetry group'
      call rval2('SYMGRPAF', ch=ch); symgaf=adjustl(ch) ! Extra Generator for adding anti ferro symmetry'
      call rval2('EWALD_AS',rr=rr,defa=[2d0]);   lat_as=rr  !'Ewald smoothing parameter
      call rval2('EWALD_TOL',rr=rr,defa=[1d-8]); lat_tol=rr !'Ewald tolerance')
      call rval2('EWALD_NKDMX',rr=rr,defa=[real(8):: 300]); lat_nkdmx=nint(rr) !'Ewald tolerance'
      mix_nsave = 8         ! nsave = # iter to save on disk
      mix_tolu = 0          ! tolu
      mix_umix = 1          ! umix (mixing parm for LDA+U)
      smalit = NULLI
      call rval2('ITER_NIT', rr=rr, defa=[real(8):: 30]); iter_maxit=nint(rr) !'maximum number of iterations in self-consistency cycle')
      call rval2('ITER_NRMIX',rr=rr,defa=[real(8):: 80]); smalit=nint(rr)     !'lmfa rseq max iter')
      call rval2('ITER_MIX', ch=ch); iter_mix= trim(adjustl(ch)) !Mixing rule for charge mixing.
      !          Anderson: A[nmix][,b=beta][,bv=betv][,n=nit][,w=w1,w2][,nam=fn][,k=nkill]
      !          Broyden:  B[nmix][,b=beta][,bv=betv][,wc=wc][,n=#][,w=w1,w2][,nam=fn][,k=nkill]
      call rval2('ITER_CONV', rr=rr,defa=[1d-4]); etol=rr   !Tolerance in energy change from prior iteration for self-consistency')
      call rval2('ITER_CONVC',rr=rr,defa=[1d-4]); qtol=rr   !Tolerance in output-input charge for self-consistency')
      call rval2('ITER_UMIX',rr=rr,defa=[.5d0]) ; mix_umix=rr !Mixing parameter for densmat in LDA+U') !2022mar9 default umix=0.5
      call rval2('ITER_TOLU',rr=rr,defa=[0d0])  ; mix_tolu=rr !Tolerance for densmat in LDA+U')
      broy  = 0
      beta  = 1d0
      wc    = -1
      wt(1:2) = 1
      wt(3) = -9 !     Flags parmxp that there are no extra elements to mix
      nmix  = -1
      if(.NOT. parmxp(trim(iter_mix),len(trim(iter_mix)),broy,nmix,wt,beta,wc,killj))& 
           call rx('MIXRHO: parse in parmxp failed')
      write(stdo,ftox)' mmmixing parameters: A/B nmix wt:',broy,nmix,ftof(wt),&
           'beta elin wc killj=',ftof(beta),ftof(wc),killj !Get wc,killj,wtinit,betainit,broyinit,nmixinit,bexist
      wtinit  =wt    !out
      betainit=beta  !out
      broyinit=broy  !out
      nmixinit=nmix  !out
      bexist=.false. 
      if(beta/=1d0) bexist=.true. !out
      call rval2('DYN_MODE',rr=rr,defa=[real(8):: 0]); lrlxr=nint(rr) ! Dynamics is only for relaxation
      !  lrlxr=  0: no relaxation, 4:relaxation(conjugate gradients), 5:relaxation(Fletcher-Powell), 6:relaxation(Broyden)
      if(lrlxr/=0) lfrce=1
      call rval2('DYN_NIT',rr=rr,  defa=[real(8):: 1]); nitrlx=nint(rr) !'maximum number of relaxation steps (statics)'//' or time steps (dynamics)')
      call rval2('DYN_HESS',rr=rr, defa=[real(8):: 1]); rdhessr= nint(rr)==1 !'Read hessian matrix')
      call rval2('DYN_XTOL',rr=rr, defa=[1d-3]);xtolr=rr !Convergence criterion in displacements XTOL>0: use length; <0: use max val; =0: do not use')
      call rval2('DYN_GTOL',rr=rr, defa=[0d0]); gtolr=rr !Convergence criterion in gradients'GTOL>0: use length;  <0: use max val;  =0: do not use')
      call rval2('DYN_STEP',rr=rr, defa=[0.015d0]); stepr=rr !Initial (and maximum) step length'
      call rval2('DYN_NKILL',rr=rr,defa=[real(8):: 0]);nkillr=nint(rr)!'Remove hessian after NKILL iter')
      !============ end of reading catok =======================      
      ! Print verbose
      i0=setprint0(30)      !initial verbose set
      if    (cmdopt2('--pr=',outs))then; read(outs,*) verbos
      elseif(cmdopt2('--pr',outs)) then; read(outs,*) verbos
      elseif(cmdopt2('-pr',outs))  then; read(outs,*) verbos
      endif
      i0 = setprint0(verbos)                   !Set initial verbos
      if( .NOT. master_mpi) i0=setprint0(-100) !iprint()=0 except master
      ! Timing: Turns CPU timing log #1:tree depth #2:CPU times as routines execute.
      !         Args may be set by command-line: --time=#1,#2
      if(cmdopt2('--time',outs) ) then 
         outs=trim(outs(2:))//' 999 999'
         read(outs,*)io_tim(1),io_tim(2)
         if(io_tim(1)==999) io_tim(1)=5
         if(io_tim(2)==999) io_tim(2)=2
         i0=1
      endif
      if(i0>=1) call tcinit(io_tim(2),io_tim(1),levelinit)
      call tcn('m_lmfinit') !after tcinit call ==============================================
      lcd4=F
      if (prgnam == 'LMF' .OR. prgnam == 'LMFGWD') lcd4=T
      if(sum(abs(socaxis-[0d0,0d0,1d0])) >1d-6 .AND. (.NOT.cmdopt0('--phispinsym'))) &
           call rx('We need --phispinsym for SO=1 and HAM_SOCAXIS/=001. Need check if you dislike --phispinsym')
      if(cmdopt0('--zmel0')) OVEPS=0d0
      if(pwmode==10) pwmode=0   !takao added. corrected Sep2011
      if(prgnam=='LMFGWD') pwmode=10+ mod(pwmode,10)
      if(iprint()>0) write(stdo,ftox) ' ===> for --jobgw, pwmode is switched to be ',pwmode
      nspecset: do 1111 j = 1, nspec
         write(stdo,'(1x)')
         write(stdo,ftox)' ... Species ',j
         !    Radial mesh parameters: determine default value of a
         i0 = NULLI
         xxx = NULLR
         call pshpr(0)
         call rmesh(z(j),rmt(j),lrel,.false.,nrmx,xxx,i0)
         call poppr
         if (xxx == .03d0) xxx = .015d0 !.025d0 jun2012 .025 to .015 as default.
         if(spec_a(j)==0d0) spec_a(j)=xxx
         i0 = 0
         call pshpr(0)
         call rmesh(z(j),rmt(j),lrel,.false.,nrmx,spec_a(j),i0)
         call poppr
         if (nr(j) == 0) nr(j) = i0
         if (rmt(j) == 0 ) lmxa(j) = -1
         lmxaj = lmxa(j)
         nlaj = 1+lmxaj
         idxdn(:,:,j) = 1
         pnuqnublock: if (nlaj /= 0) then
            !! ==== Reset default P,Q in absence of explicit specification ====
            !     ! -- takao jun2012. qnu is set by default p. --
            !     ! This looks too complicated. Fix this in future.
            !     ! In anyway, we expect pnusp and qnu are correctly returned (qnu does not care value of given P).
            !     print *,'qnuin ',sum(abs(qnu(:,:,j))),qnu(:,:,j)
            !     ! set default pnusp. See the following section 'correct qnu'
            !     ! isp=1 means charge. isp=2 means mmom
            if(allocated(pnuspdefault)) deallocate(pnuspdefault,qnudefault,qnudummy)
            allocate(pnuspdefault(n0,nsp),qnudefault(n0,nsp),qnudummy(n0,nsp))
            pnuspdefault=0d0
            qnudefault=0d0
            qnudummy=0d0
            iqnu=1
            if(sum(abs(qnu(:,1,j)))<1d-8) iqnu=0 !check initial Q is given or not.
            call defpq(z(j),lmxaj,1,pnuspdefault,qnudefault)! qnu is given here for default pnusp.
            call defpq(z(j),lmxaj,1,pnusp(1,1,j),qnudummy)  ! set pnusp. qnu is kept (but not used here).
            if(iqnu==0) qnu(:,1,j)=qnudefault(:,1)
            if (nsp == 2) call dcopy(n0,pnusp(1,1,j),1,pnusp(1,2,j),1)
            if(nsp==2) pzsp(1:n0,2,j) = pzsp(1:n0,1,j) !takao
            !! lmxb corrected by pzsp
            nnx=0 !nout
            do i=n0,1,-1
               if(pzsp(i,1,j)>0d0) then
                  nnx=i
                  lmxb(j)=max(lmxb(j),nnx-1)
                  exit
               endif
            enddo
            if (nnx>0) then
               if (maxval(pzsp(1:nnx,1,j))>0) then
                  lpzi = 1 !max(lpzi,1) !,2
                  lpz(j)=1
                  if ( sum(floor(pzsp(1:nlaj,1,j)/10))>0 ) then
                     lpzex(j)=1
                  endif
               endif
            endif
            if(maxval(pzsp(1:nnx,1,j))>10d0) lpztail= .TRUE. ! PZ +10 mode exist or not.
            readpnublock:block ! P = PrincipleQnum - 0.5*atan(dphidr/phi)/pi
              integer:: ifipnu,lr,iz,nspx,lrmx,isp,ispx
              real(8):: pnur,pzav(n0),pnav(n0)
              character(8):: charext
              if (prgnam /= 'LMFA'.and.ReadPnu) then
                 open(newunit=ifipnu,file='atmpnu.'//trim(charext(j))//'.'//trim(sname))
                 write(stdo,*)'READP=T: read pnu from atmpnu.*'
                 do
                    read(ifipnu,*,end=1015) pnur,iz,lr,isp
                    if(iz==1) pzsp (lr+1,isp,j)= pnur ! +10d0 caused probelm for 3P of Fe.
                    if(iz==0) pnusp(lr+1,isp,j)= pnur
                    lrmx=lr
                    nspx=isp
                 enddo
1015             continue
                 pzav(1:lrmx+1)=sum(pzsp(1:lrmx+1,1:nspx,j), dim=2)/nspx !spin averaged
                 pnav(1:lrmx+1)=sum(pnusp(1:lrmx+1,1:nspx,j),dim=2)/nspx
                 do l=1,lrmx+1
                    if(pzav(l)>pnav(l)) pzav(l)=floor(pzav(l))+.5d0 ! push up p =floor(p)+0.5 if we have lower orbital
                    if(pzav(l)>1d-8.and.pnav(l)>pzav(l)) pnav(l)=floor(pnav(l))+.5d0
                 enddo
                 do ispx=1,nspx
                    pzsp(1:lrmx+1, ispx,j) = pzav(1:lrmx+1)
                    pnusp(1:lrmx+1,ispx,j)=  pnav(1:lrmx+1)
                 enddo
                 close(ifipnu)
              endif
            endblock readpnublock
            ! our four cases are
            !     P=Pdefault      ! qnu
            !     Pdefault < P    ! Pdefault is filled as core
            !     Pz < P=Pdefault ! qnu + 2*(2l+1)
            !     Pz=Pdefault < P ! qnu
            if(iqnu==0) then
               do lx=0,lmxaj     !correct valence number of electrons.
                  if(pzsp(lx+1,1,j)<1d-8) then ! PZSP(local orbital) not exist
                     if( int(pnuspdefault(lx+1,1)) < int(pnusp(lx+1,1,j)) ) then
                        qnu(lx+1,1,j)= 0d0 ! pnuspdefault is filled and no q for pnusp. (core hole case or so)
                     endif
                  else           !PZ exist   !     print *,'qnu=',lx,qnu(lx+1,1,j)
                     if( mod(int(pzsp(lx+1,1,j)),10)<int(pnuspdefault(lx+1,1)) ) then
                        qnu(lx+1,1,j)= qnu(lx+1,1,j)+ 2d0*(2d0*lx+1d0)
                     endif
                  endif
               enddo
            endif
         endif pnuqnublock    
         rg(j)= 0.25d0*rmt(j)
         rfoca(j)= 0.4d0*rmt(j)
         if(procid==master .AND. sum(abs(idu(:,j)))/=0) then
            inquire(file='sigm.'//trim(sname),exist=sexist)
            if(sexist) then
               do lxx=0+1,3+1
                  if(idu(lxx,j)>10) then
                     write(stdo,"(a,2i4)")'For IDU>10 with sigm.*, we set UH=JH=0 for l,ibas=',lxx,j
                     uh(lxx,j) = 0d0
                     jh(lxx,j) = 0d0
                  endif
               enddo
            endif
            do lxx=0+1,3+1
               idu(lxx,j) = mod(idu(lxx,j),10)
            enddo
         endif
         call mpibc1_int(idu(:,j),4,'m_lmfinit_idu')
         call mpibc1_real(uh(:,j),4,'m_lmfinit_uh') !bug fix at oct26 2021 (it was _int-->not passed )
         call mpibc1_real(jh(:,j),4,'m_lmfinit_jh')
         nkaphh(j) = nkapii(j) + lpzex(j)
1111  enddo nspecset
      lmxax = maxval(lmxa) !Maximum L-cutoff
      nlmax = (lmxax+1)**2
      nkaph = nkapi + lpzi     !-1
      mxorb= nkaph*nlmax
    endblock Stage1GetCatok
    Stage2SetModuleParameters: block
      integer:: isw,iprint
      logical:: cmdopt0
      maxit=iter_maxit
      nl = max(lmxbx,lmxax)+1 !max l-base l-aug +1
      nlmax=nl**2
      if (dalat == NULLR) dalat=0d0
      lat_alat=alat+dalat
      lat_avw=avw
      lat_nkqmx=lat_nkdmx
      lat_platin=plat
      lat_tolft=tolft
      !! set cg coefficients for lmf --- Choose dimensions for arrays
      lmxcg=8
      lmxcy=12
      !      if (lmxcg .le. 6) then
      !        lnjcg = 6500
      !        lnxcg = 1300
      !      else if (lmxcg .le. 8) then
      lnjcg = 22700
      lnxcg = 3400
      !      else if (lmxcg .le. 10) then
      !        lnjcg = 62200
      !        lnxcg = 7400
      !      else
      !        call rxi('setcg: cannot handle lmxcg=',lmxcg)
      !      endif
      nlm=(lmxcy+1)**2
      allocate(rv_a_ocy(abs(nlm)))
      allocate(rv_a_ocg(abs(lnjcg)))
      allocate(iv_a_ojcg(abs(lnjcg)))
      allocate(iv_a_oidxcg(abs(lnxcg)))
      call sylmnc ( rv_a_ocy , lmxcy )
      call scg ( lmxcg , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg ) !set CG coefficients for lmf part.
      ham_nkaph=nkaph
      if (procid==master) then
         inquire(file='sigm.'//trim(sname),exist=sexist)
         if (lrsigx/=0 .AND. ( .NOT. sexist) ) then
            write(stdo,*)' bndfp (warning): no sigm file found ... LDA calculation only'
            lrsigx = 0
         endif
      endif
      ham_lsig=lrsigx
      call mpibc1_int(ham_lsig,1,'bndfp_ham_lsig')
      ham_scaledsigma=scaledsigma
      ham_pwmode=pwmode
      ham_oveps=oveps
      !! idxdn= 1 or 3
      !r     1  Orbital is included as "active" orbitals, which means
      !r           they are included in the hamiltonian and diagonalized
      !r     3   Neglected.
      ! xxxxxxxxx
      !r     10  Orbital is a local orbital whose value and slope are
      !r         constructed to be zero at the MT boundary.
      !r         It is included in the basis.
      !r     11  Orbital is a local orbital with a smooth Hankel tail
      !r         and it is included in the basis.
      do j=1,nspec
         do  ik = 1, nkap0
            if (ik <= nkapi) then
               do lp1 = 1, lmxb(j)+1
                  if (ik==1 .AND. rsmh1(lp1,j)<=0) idxdn(lp1,ik,j)=3
                  if (ik==2 .AND. rsmh2(lp1,j)<=0) idxdn(lp1,ik,j)=3
               enddo
            endif
            if (ik == nkaph .AND. sum(lpz)>0) then
               idxdn(:,ik,j)=3  !call ivset(idxdn(1,ik,j),1,n0,4)
               do  lp1  = 1, lmxb(j)+1
                  if (pzsp(lp1,1,j) /=  0) then
                     if(pzsp(lp1,1,j)>=10) idxdn(lp1,ik,j)=1 !11
                     if(pzsp(lp1,1,j)>  0) idxdn(lp1,ik,j)=1 !10
                  endif
               enddo
            endif
            if (ik > nkaph) idxdn(:,ik,j)=3 !call ivset(idxdn(1,ik,j),1,n0,4)
            idxdn(lmxb(j)+ 2:,ik,j)=3
         enddo
      enddo
      allocate(v_sspec(nspec))
      do j=1,nspec !additional data supplied from rdovfa.f90 and iors.f90
         v_sspec(j)%z=z(j)
         v_sspec(j)%a=spec_a(j)
         v_sspec(j)%nr=nr(j)
         v_sspec(j)%kmxt=kmxt(j)
         v_sspec(j)%lfoca=lfoca(j) !lfoca=1,usually (frozen core)
         v_sspec(j)%rsmv= rmt(j)*.5d0 !rsmv(j)
         v_sspec(j)%lmxa=lmxa(j) !lmx for augmentation
         v_sspec(j)%lmxb=lmxb(j) !lmx for base
         v_sspec(j)%lmxl=lmxl(j) !lmx for rho and density
         v_sspec(j)%rfoca=rfoca(j)
         v_sspec(j)%rg=rg(j)
         v_sspec(j)%rmt=rmt(j)
      enddo
      do j=1,nbas
         is=ispec(j) 
         pnuall(:,1:nsp,j) = pnusp(1:n0,1:nsp,is)
         pnzall(:,1:nsp,j) = pzsp(1:n0,1:nsp,is)
         if(procid==master) then
         do isp=1,nsp
         write(6,ftox)'pnuall: j isp pnu=',j,isp,ftof(pnuall(1:lmxa(is)+1,isp,j),6)
         write(6,ftox)'pnzall: j isp  pz=',j,isp,ftof(pnzall(1:lmxa(is)+1,isp,j),6)
         enddo
         endif
      enddo
      !! ... Suppress symmetry operations for special circumstances
      ! addinv=T only when we expect Hamiltonian is real even with keeping spin direction.
      ! addinv=T means psi* is eigenfunction.(Time-Reversal with keeping spin).
      ! 2022-jan-20 new setting of addinv (addinv =.not.ctrl_noinv)
      !Add inversion to get !When we have TR, psi_-k(r) = (psi_k(r))^* (when we have SO/=1).
      !                      density |psi_-k(r)|^2 = |psi_k^*(r)|^2
      if((lrlxr>=1 .AND. lrlxr<=3) .OR. &
           cmdopt0('--cls') .OR. cmdopt0('--nosym') .OR. cmdopt0('--pdos')) then
         symg = 'e'
         addinv = .false. !add inversion 
      elseif(lso == 0) then
         addinv=.true. !add inversion means
      else
         addinv=.false. 
      endif
      sstrnsymg=trim(symg)
      nspc = 1
      if( lso==1 ) nspc = 2
      orbital: block
        !o   norb  :number of orbital types for ib; see Remarks
        !o   ltab  :table of l-quantum numbers for each type
        !o   ktab  :table of energy index for each type
        !o   offl  :offl(norb) offset in h to this block of orbitals
        !o   ndim  :dimension of hamiltonian for this site
        !r Remarks
        !r   Each orbital type is label by a 'l' and a 'k' index
        !r   Each orbital corresponds to a unique radial wave function at the
        !r   site where the orbit is centered.  There can be multiple 'k'
        !r   indices (radial wave function shapes) for a particular l.
        integer:: ib,l,lmr,ia, nnrlx,lmri,ik,nnrl,nnrli,li, iprmb(nbas * nl**2 * maxp )
        iprmb=-1
        nlmto = 0
        do 110 ib = 1,nbas
           is = ispec(ib) 
           iposn = mxorb*(ib-1)
           do 1121 ik = 1, nkaph
              do  l = 0, nl-1
                 do  m = -l, l
                    iposn = iposn+1
                    if(idxdn(l+1,ik,is)==1) then
                       nlmto = nlmto+1
                       iprmb(iposn) = nlmto
                    endif
                 enddo
              enddo
1121       enddo
110     enddo
        nnrlx=0
        do ib = 1,nbas
           nnrli = 0
           lmri = nl*nl*nkaph*(ib-1)
           do  ik = 1, nkaph
              do   li = 0, nl-1
                 lmri = lmri + 2*li+1
                 if (iprmb(lmri) > nlmto) cycle
                 nnrli = nnrli + 2*li+1
              enddo
           enddo
           nnrlx = max(nnrlx,nnrli)
        enddo
        nnrl = nnrlx
        ndimx=nnrl
        allocate(ltabx(n00,nbas),ktabx(n00,nbas),offlx(n00,nbas),ndimxx(nbas),norbx(nbas))
        norbx=0
        ndimxx=0
        do ib=1,nbas
           lmr = nl*nl*nkaph*(ib-1)
           do  ik = 1, nkaph
              do  l = 0, nl-1
                 offlx(norbx(ib)+1,ib) = -1
                 if (iprmb(lmr+1) >0 .AND. iprmb(lmr+1) <= nlmto) then
                    offlx(norbx(ib)+1,ib) = iprmb(lmr+1) - 1
                    norbx(ib) = norbx(ib) +1
                    ndimxx(ib) = ndimxx(ib) + 2*l+1
                    ltabx(norbx(ib),ib) = l
                    ktabx(norbx(ib),ib) = ik
                 endif
                 lmr = lmr + 2*l+1
              enddo
           enddo
           if (norbx(ib) > n00) call rx('orbl: norb> n00')
        enddo
      endblock orbital
      nspx  = nsp
      if(lso/=0) nspx = 1
      nvi= sum([( (lmxa(ispec(ib))+1)**2,ib=1,nbas )])
      nvl= sum([( (lmxl(ispec(ib))+1)**2,ib=1,nbas )])
      pot_nlma=nvi
      pot_nlml=nvl
      allocate(jnlml(nbas),jma(nbas))
      jnlml(1)=1 ! (ilm,ib) index
      jma(1)=1
      do  i = 1, nbas-1
         jnlml(i+1)= (lmxl(ispec(i))+1)**2 +jnlml(i)
         jma(i+1)= (lmxa(ispec(i))+1)**2 +jma(i)
      enddo
      !     Make nat = number of real atoms as nbas - # sites w/ floating orbitals
      if (procid == master) then
         nat = nbas
         do  i = 1, nbas
            j=ispec(i) 
            l=v_sspec(j)%lmxa
            if (l == -1) nat = nat-1
         enddo
      endif
      call mpibc1_int(nat,1,'m_lmfinit_nat')
      allocate(wowk(nbas))
      wowk=0
      call pshpr(0)
      call suldau(nbas,nlibu,k,wowk)!Count LDA+U blocks (printout only)
      ham_nlibu=nlibu
      call poppr
      deallocate(wowk,amom, stni,rg,rfoca,idxdn, rmt,  lfoca,lmxl, spec_a,nr) !,rsmv)
      !! --- takao embed contents in susite here. This is only for lmf and lmfgw.
      allocate(iv_a_oips(nbas),source=[(ispec(ib), ib=1,nbas)])
      seref= sum([(eref(ispec(ib)),ib=1,nbas)])
      ham_seref= seref
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if( cmdopt0('--quit=show') ) call rx0(trim(prgnam)//' --quit=show')
    endblock Stage2SetModuleParameters
    Stage3InitialSetting :block 
      logical:: cmdopt0
      integer:: iprint
      if (cmdopt0('--etot')) then
         lfrce=0
         maxit=1
      endif
      if(lhf) maxit= 1
      if (lrel /= 0) then
         cc = 274.074d0 !srel
      else 
         cc = 1d10      !nrel 
      endif
      call setcc(lrel) !lrel/=0 means scalar relativistiv c=274.074d0 in a.u.
      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') !fullmesh stop just after do 2010 iq loop.
      lrout = 1 ! Whether to evaluate output density and/or KS energy
      leks = 1
      if (bz_nevmx == -1) then !bz_nevmx ! Maximum number of eigenvalues
         lrout = 0
         leks = 0
      endif
      if(cmdopt0('--band') .OR. fullmesh) then !Switch to plot bands at specified qp
         lfrce = 0
         plbnd = 1
         lrout = 0
      else
         plbnd=0
      endif
      if(lrout == 0 ) maxit = 1
      if(bz_lmet/=0 .OR. bz_tetrahedron ) ldos=1
      if(ldos==0) bz_ndos = 1
      if(lrout == 0 .AND. lfrce /= 0) then
         write(stdo,"(a)") 'Error: output density required when forces sought.\n'// &
              '      To make output density turn off HF=t and/or NEVMX<0'
         call Rx('incompatible input. see the end of console output')
      endif
      if(lrout == 0 .AND. cmdopt0('--etot')) then
         write(stdo,"(a)") 'Error: output density required with --etot switch.\n'// &
              '      To make output density turn off HF=t and/or NEVMX<0'
         call Rx('incompatible input. see the end of console output')
      endif
      !! LDA+U block (it was suldau.F)
      nlibu = 0
      lmaxu = 0
      allocate(lldau(nbas))
      lldau = 0
      if(master_mpi) write(stdo,*)
      do  ib = 1, nbas
         is = ispec(ib)
         do  lx = 0, min(v_sspec(is)%lmxa,3)
            if (idu(lx+1,is) /= 0) then
               if (lldau(ib) ==0) lldau(ib) = nlibu+1
               nlibu = nlibu+1
               lmaxu = max(lmaxu,lx)
               !if(master_mpi)write(stdo,ftox)'lda+u block ibas lx=', ib,lx, &
               !'JH=',ftof(v_sspec(is)%jh(lx+1),3),'UH',ftof(v_sspec(is)%uh(lx+1),3)
            endif
         enddo
      enddo       !! aug2012 we now fix lcplxp=1 (complex ppi integral)
      ! lhh, nkapii, nkaphh (nkaphh = nkapii(1 or 2) +1) if extented local orbital exist)
      allocate(lhh(nkap0,nspec))
      lhh=-1
      do i=1,nspec
         lmxbj = lmxb(i)
         call getiout(rsmh1(1,i), lmxbj+1,lhh(1,i))
         if(nkapii(i)==2) call getiout(rsmh2(1,i),lmxbj+1,lhh(2,i))
         if(lpz(i)==1 )   call getiout(pzsp(1,1,i),lmxbj+1,lhh(nkaph,i))!lh for lo
      enddo
      if(master_mpi) then
         write(stdo,*)
         write(stdo,"('mmm === MTO setting ===')")
         do i=1,nspec
            lmxbj = lmxb(i)
            write(stdo,"('mmm ispec lmxb lpz nkapii nkaphh=',10i5)")i,lmxb(i),lpz(i),nkapii(i),nkaphh(i)
            write(stdo,"('mmm rsmh1 ',i4,100f6.2)")i, rsmh1(1:lhh(1,i)+1,i)
            write(stdo,"('mmm   eh1 ',i4,100f6.2)")i,   eh1(1:lhh(1,i)+1,i)
            if(nkapii(i)==2) write(stdo,"('mmm rsmh2 ',i4,100f6.2)")i, rsmh2(1:lhh(2,i)+1,i)
            if(nkapii(i)==2) write(stdo,"('mmm  eh2  ',i4,100f6.2)")i,   eh2(1:lhh(2,i)+1,i)
            if(lpz(i)==1 ) write(stdo,"('mmm pz    ',i4,100f6.2)")i,    pzsp(1:lhh(nkaph,i)+1,1,i)
            write(stdo,"('mmm lh    ',i4,100i3)")  lhh(1:nkaph,i)
         enddo
      endif
      relaxmodesetting: block ! Atomic position Relaxation setup (DYN mode)
        integer:: i,j,k,iprint !,ifrlx(3)
        logical:: force,mdxx
        ! nitrlx = num of iteration cycle for atomic relaxiation (outer loop)
        if(lrlxr==0)nitrlx=0 !no relaxiation. Only sc calculation for given atomic position.
        if(lrlxr>0 .AND. lrlxr<4 ) call rx('lmf not set up for MD yet')
        if(nitrlx>0) then !nitrlx >0 is for atomic position relaxiation lmf.f90
           allocate(indrx_iv(2,3*nbas))
           force = lfrce>0
           if ( .NOT. force .OR. lrlxr== 0) goto 9299
           mdxx =  lrlxr <= 3
           j = 0
           if (mdxx) then
              xyzfrz = .false.
              goto 9299
           elseif (force) then
              do  i = 1, nbas
                 do  k = 1, 3
                    if (ifrlx(k,i) == 1) then
                       j = j + 1
                       indrx_iv(1,j) = k
                       indrx_iv(2,j) = i
                       xyzfrz(k) = .false.
                    endif
                 enddo
              enddo
           endif
           natrlx = j
           if (natrlx == 0) goto 9299
           pdim = 0
           if ( .NOT. mdxx) then
              if (lrlxr == 4) pdim = natrlx*7
              if (lrlxr == 5) pdim = natrlx*(7+natrlx)
              if (lrlxr == 6) pdim = natrlx*(12+2*natrlx)
           endif
           if (iprint() >= 30) then
              if (lrlxr == 4) then
                 write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (conjugate gradients) ..'
              elseif (lrlxr == 5) then
                 write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (Fletcher-Powell) ..'
              else
                 write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (Broyden) ..'
              endif
              write(stdo,ftox)'   relaxing',natrlx,'variables,',nitrlx,'iterations'
              write(stdo,ftox)'   x-tol g-tol step=', xtolr,gtolr,stepr
           endif
9299       continue
        endif
      endblock relaxmodesetting
    endblock Stage3InitialSetting
    call tcx('m_lmfinit')
  end subroutine m_lmfinit_init
  subroutine getiout(a,iin,iout) !a(1:iout) can be nonzero.
    integer:: iin,iout,i
    real(8):: a(iin)
    do i=iin,1,-1
       if(a(i)>0) then
          iout=i-1 !angular mom
          exit
       endif
    enddo
  end subroutine getiout
  subroutine fill3in(nin,res)!- fills res(2) or res(2:3) if res(1) or res(1:2) are given
    integer :: nin,res(3)
    if (nin==2) then
       res(3) = res(2)
    elseif (nin==1) then
       res(2:3) = res(1)
    endif
  end subroutine fill3in
end module m_lmfinit

! mmmm old doc mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
! Hereafter are list of variables. Old document, but it maybe a help.
! TK tried to remove obsolate things 
! def      Uncertainty in Fermi level
! dosw     Energy window over which DOS accumulated
! fsmom    fixed-spin moment (fixed spin moment method)
! ef       Fermi level
! ntet     number of tetrahedra
! efmax    eigenvalues above efmax are not needed
! lmet     0 insulator   
!          3 metal
!r          (numbers, not powers of 2)
!r lshft    shift mesh in each of three directions
!r bz_n     Polynomial order for Methfessel-Paxton sampling
!r ndos     No. DOS points (sampling integration, and lmdos)
!r nevmx    eigenvectors above nevmx are not needed
!r nkabc    number of divisions in qp mesh
!r bz_w        Line broadening for sampling integration
!
!   lcd4  represent full potential density on a uniform mesh
!   ldos    1 make dos
!r   lfrce   How forces are calculated
!r           0 do not calculate forces
!r           1 shift in FA density
!r   noinv    1 do not add inversion !noinv BZ_NOINV
!r   lrel    1 scalar relativistic
!r           2 fully relativistic
!r   lxcf    parameter defining XC functional
!r           1s digit:
!r           1 for Ceperly-Alder
!r           2 for Barth-Hedin (ASW fit)
!r           103 for PBE-GGA
!r   maxit   max. no.  iterations in self-consistency cycle
!r   nclass  size of class
!r   nl      1 + maximum lmxa
!r   nbas   number of sites
!r   nspec   number of species
!r   nspin   number of spins
!r   omax1   sphere overlap constraints, type 1
!r   omax2   sphere overlap constraints, type 2
!r   rmaxes  upper limit to ES radius when finding new empty spheres
!r   rmines  lower limit to ES radius when finding new empty spheres
!r   sclwsr  scale wsr until reaching this fractional vol
!r           10s digit used for assymetric treatment of ES:
!r             0 ES and other sites are treated symmetrically
!r             1 all sites with z>0 are resized first; then
!r               all sites are resized.
!r             2 all sites with z>0 are resized first; then
!r               the ES sites only are resized
!r   tol     1 q- tolerance for self-consistency
!r           2 e- tolerance for self-consistency
!r   wsrmax  constraint on size of largest WS sphere
!r   zbak    background charge and MT correction parameters
!r                     zbak = uniform background charge included
!r                            in electrostatics but not in xc pot.
!r  ehf    Harris-Foulkes energy
!r  ehk    Hohnberg-Kohn energy
!r  nkaph  number of repetitions of one l-quant.n. in basis
!r  nlibu  Number of LDA+U blocks
!
!r  pwemax High Energy cutoff for PW part of basis
!r  pwmode Controls PW part of basis
!r         0 => no PW part of basis
!r         1 => include PWs in basis
!r         2 => include only PWs in basis   +10 means |G| truncation.
!
!r  alat    lattice parameter, in a.u.
!r  avw     average MT radius
!r  awald   Ewald smoothing parameter
!r  nabc    no. divisions for F.T. mesh
!r  ng      no. G vectors
!r  nkd     no. direct latt. vecs. for Ewald sum
!r  nkdmx   dimensioning for arrays holding latt. vecs
!r  nkq     no. reciprocal latt. vecs. for Ewald sum
!r  nkqmx   dimensioning for arrays holding latt. vecs
!r  npgrp   Number of point symmetry group operations
!r  nsgrp   Number of space symmetry group operations
!r  plat..  lattice vectors, units of alat
!r  qlat..  reciprocal lattice vectors, units 2pi/a
!r  rpad..  truncate Ewald to rpad*rmax when lattice vector
!r          list has to be padded in order to include at
!r          least one lattice vector
!r  tol     Ewald tolerance
!r  tolft   FT mesh tolerance
!r  vol     cell volume
!
!r   b       mixing beta
!r   bv      extra potential mixing
!r   fn      mixing file name
!r   kill    kill the mixing file after k iterations
!r   lxpot   decouple potential and the charge
!r           1: mix vin and v(qmix)
!r           2: mix vin and v(qout)
!r   mmix    maximum number to mix
!r   mode    1 Anderson 2 Broyden
!r   n       Number of iterations for this species
!r   nsave   number of iterations to save on disk
!r   r..     expression for rmscst
!r   rms1    1st rms error
!r   rms2    2nd rms error
!r   size    total size of this struc
!r   tj..    Anderson t's
!r   tolu    tolerance for LDA+U
!r   umix    mixing parameter for LDA+U
!r   w..     Linear mixing weights
!r   wc      Broyden weight



