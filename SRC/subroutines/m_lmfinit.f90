!> Ititial data for lmf-MPIK lmchk lmfa read from ctrl file
!> pos can be from AtomPos (see lmfp.f90)
module m_lmfinit ! 'call m_lmfinit_init' sets all initial data from ctrl are processed and stored in m_lmfinit_init.
  use m_ftox !for write(*,ftox) ftof(values)
  use m_ext,only :sname        ! sname contains extension. foobar of ctrl.foobar
  use m_MPItk,only: master_mpi
  use m_lgunit,only: stdo,stdl
  use m_fatom,only:   sspec !allocation only in m_lmfinit: free atom density (detremined by lmfa) WARN: Not protected.
  use m_density,only: pnuall,pnzall !NOTE: These are set here! log-derivative of radial functions. m_denisty is NOT protected.

  use m_nvfortran,only: findloc
  
  implicit none 
  public:: m_lmfinit_init,icgi,icge
  integer,public,parameter:: noutmx=48,NULLI=-99999,nkap0=3,mxspec=256,lstrn=1000,n0=10,nppn=2,nrmx=1501,nlmx=64,n00=n0*nkap0,k0=3
  real(8),public,parameter:: fpi=16d0*datan(1d0), y0=1d0/dsqrt(fpi), pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi),pi4=fpi,&
       NULLR =-99999, fs = 20.67098d0, degK = 6.3333d-6 ! defaults for MD
  logical,public,parameter:: T=.true., F=.false.
  integer,public,protected:: lat_nkqmx,lat_nkdmx, lxcf, smalit,lstonr(3)=0,nbas=NULLI,nspec,&
       nspc,master=0,nspx, maxit,gga,ftmesh(3),nmto=0,lrsigx=0,nsp=1,lrel=1,lso=0,&
       lmxbx=-1,lmxax,bz_lshft(3)=0, bz_lmet,bz_n,bz_lmull,bz_fsmommethod,str_mxnbr,&
       iter_maxit=1, mix_nsave, pwmode,ncutovl ,ndimx,natrlx, leks,lrout,plbnd, pot_nlma, pot_nlml,ham_nspx, nlmto,& !total number of MTOs 
       lrlxr,nkillr,nitrlx, broyinit,nmixinit,killj ,&
       ham_pwmode,ham_nlibu, nlmax,lfrce,bz_nevmx,ham_nbf,ham_lsig,bz_nabcin(3)=NULLI, bz_ndos,ldos, lmaxu,nlibu,nlmxlx
  logical,public,protected :: ham_frzwf,ham_ewald, lhf,lcd4,bz_tetrahedron, addinv,&
       readpnu,v0fix,pnufix,bexist,rdhessr, lpztail=.false., readpnuskipf, afsym
  real(8),public,protected:: pmin(n0)=0d0,pmax(n0)=0d0,tolft,scaledsigma, ham_oveps,ham_scaledsigma, cc,&!speed of light
       dlat,alat=NULLR,dalat=NULLR,vol,avw,vmtz(mxspec)=-.5d0,&
       bz_efmax,bz_zval,bz_fsmom,bz_semsh(10),zbak,bz_range=5d0,bz_dosmax,&
       lat_as,lat_tol,lat_rpad=0d0, str_rmax=nullr, etol,qtol,mix_tolu,mix_umix, pwemax,oveps,& !,pwemin=0d0
       ham_seref, bz_w,lat_platin(3,3),lat_alat,lat_avw,lat_tolft,lat_gmaxin,socaxis(3), xtolr,gtolr,stepr, wtinit(2),wc,betainit 
  character(lstrn),public,protected:: sstrnsymg, symg=' ',   symgaf=' '!for Antiferro
  character(128),public,protected :: iter_mix=' ' !mix
  character(8),public,protected:: alabl
  character(8),public,allocatable,protected:: slabl(:)
  logical,public,allocatable,protected:: cstrmx(:),frzwfa(:)
  integer,public,allocatable,protected:: lmxb(:),lmxa(:),idmod(:,:),idu(:,:),kmxt(:),lfoca(:),lmxl(:),nr(:),nmcore(:),&
       nkapii(:),nkaphh(:),ispec(:),ifrlx(:,:),ndelta(:), iantiferro(:), iv_a_oips (:),  lpz(:),lpzex(:),lhh(:,:) ,jnlml(:)
  real(8),public,allocatable,protected:: rsmh1(:,:),rsmh2(:,:),eh1(:,:),eh2(:,:),rs3(:),alpha(:,:),uh(:,:),jh(:,:),eh3(:),&
       qpol(:,:),stni(:),pnusp(:,:,:),qnu(:,:,:),pnuspdefault(:,:),qnudefault(:,:),qnudummy(:,:),&
       coreq(:,:),rg(:),rsma(:),rfoca(:), rmt(:),pzsp(:,:,:), amom(:,:),spec_a(:),z(:),eref(:),rsmv(:)
  character*(8),public,allocatable,protected:: coreh(:)
  real(8),public,allocatable,protected:: pos(:,:), delta(:,:),mpole(:),dpole(:,:)
  !WARNINIG: pos is initial position read from ctrl. We use lattic.f90 (search pos in bndfp.f90 and setpos in lmfp.f90 as well.)
  real(8),public,allocatable,protected:: rv_a_ocg (:), rv_a_ocy(:)  !ClebshGordon coefficient (GW part use clebsh_t)
  integer,public,allocatable,protected:: iv_a_oidxcg(:),iv_a_ojcg(:)!ClebshGordon coefficient (GW part use clebsh_t)
  !functions icgi,icge are for ClebshGordon coefficients. To understand how to use CG coefficients, See samples in m_smhankel.f90 and so on.
  integer,public,allocatable,protected:: lldau(:), indrx_iv(:,:) 
  integer,public,allocatable,target:: ltabx(:,:),ktabx(:,:),offlx(:,:),ndimxx(:),norbx(:),blksx(:,:),ntabx(:,:)
  integer,public,protected :: lxx,kxx,norbmto,norbxx !oribtal index
  integer,public,allocatable,protected:: ib_table(:),k_table(:),l_table(:),ltab(:),ktab(:),offl(:), offlrev(:,:,:),ibastab(:)
  private
contains
  subroutine m_lmfinit_init(prgnam,commin) ! All the initial data are set in module variables from ctrlp.*
    use m_gtv2,only: gtv2_setrcd,rval2
    use m_defpq,only:defpq
    ! Inputs
    !   file  : read ctrl.sname
    !   prgnam: name of main program
    ! Outputs
    !    All the module variables. Only several components of v_sspec are added by iors/rdovfa (readining atomic or previous results).
    !MEMO:2023-sep
    ! Note our block coding: Search ConvertCtrl2CtrlpByPython ReadCtrlp Stage1 Stage2 Stage3.
    !   At ConvertCtrl2CtrlpByPython, we convert ctrl.foobar to ctrlp.foobar by invoking a python script.
    !   BZ_  : Brillouin Zone related
    !   HAM_ :  Hamiltonian related
    !   SITE_: site information
    !   slabl : vector of species labels (spec < class < nbas)
    !   avw   : the average Wigner-Seitz radius
    !   lrel  :specifies type of Schrodinger equation
    !         :0 nonrelativistic Schrodinger equation
    !         :1 scalar relativistic Schrodinger equation
    !   lxcf  :specifies type of XC potential. 
    !         :1 for Ceperly-Alder
    !         :2 for Barth-Hedin (ASW fit)
    !         :103 for PBE
    !   
    !         z=     z(j)      !nucleus
    !         a=     spec_a(j) !a for radial mesh
    !         nr=    nr(j)     !nr for radial mesh
    !         rmt=   rmt(j)    ! MT radius
    !         rsmv=  rmt(j)*.5d0 !rsmv(j) :smoothing radius for P_kl expansion
    !         kmxt=  kmxt(j) !kmax the max number radial funciton index of P_kl(r)
    !         lmxa=  lmxa(j) !lmx for augmentation
    !         lmxb=  lmxb(j) !lmx for basis
    !         lmxl=  lmxl(j) !lmx for rho and density
    !         rfoca= rfoca(j) !smoothing radius for frozen core overlap approx
    !         lfoca= lfoca(j) !lfoca=1,usually (frozen core mode)
    !         rg=    rg(j)   !rsm for gaussians to fix multipole moments
    
    !  MTO is specified by (n,l,m). (n=1,2,3. n=1:EH1, n=2:EH2, n=3:PZ)
    !   nbas  :number of atoms in the basis
    !   nkaphh :The maximum number of radial functions centered at particular R and l channel used in the lmto basis. 
    !   nsp   :number of spins
    !   nspc  :2 if two-spins are coupled (spin-off diagonal included). Thus isp runs isp=1,nsp/nspc, where nspx=nsp/nspc
    !   nspec :number of species
    !   stdl  :standard log file
    !   stdo  :standard output file
    ! DYN_ molecular dynamics section DYN (only relaxiation, 2022-6-22)
    !   lrlxr: 0 no relaxation or dynamics, 4 relax with conjugate gradients, 5 relax with variable metric, 6 relax with Broyden
    !   rdhessr: T read hessian matrix, xtolr: relaxation x-tolerance, gtolr: relaxation g-tolerance, stepr: step length
    !   nkillr: Remove hessian after this many steps
    implicit none
    include "mpif.h" 
    character,intent(in)::  prgnam*(*)
    integer,parameter:: recln=512
    character:: strn*(recln),strn2*(recln)
    character(:),allocatable:: recrd(:)
    integer,parameter:: maxp=3
    character fileid*64
    character(256)::  a,outs,sss,ch
    character(128) :: nm
    logical :: cmdopt0,  debug,sexist=.false., ipr10,fullmesh,lzz, logarr(100)
    integer :: i,is,iprint, iprt,isw,ifi,j,k,l,lfrzw,lrs,k1,k2,lmxbj,lmxaj,nlbj,&
         ibas,ierr,lc, iqnu=0, ifzbak,nn1,nn2,nnx,nlaj,isp,&
         inumaf,iin,iout,ik,iprior,ibp1,indx,iposn,m,nvi,nvl,nn1xx,nn2xx, nnn,ib,&
         lmxcg,lmxcy,lnjcg,lnxcg,nlm,nout,nn,i0,ivec(10),iosite,io_tim(2),verbos,&
         lp1,ii,sw,it,levelinit=0, lx,lxxx, reclnr,nrecs,nrecs2 ,nkapi!,mxorb
    real(8):: pnuspx(20) ,temp33(9),seref, xxx, avwsr, d2,plat(3,3),rydberg,rr, vsn,vers,xv(2*n0),xvv(3)
    real(8),allocatable ::rv(:)
    character*(8),allocatable::clabl(:)
    integer,allocatable:: idxdn(:,:,:) 
    real(8),allocatable:: pnuspc(:,:,:),qnuc(:,:,:,:),pp(:,:,:,:),ves(:),zc(:) !    debug = cmdopt0('--debug')
    integer,optional:: commin
    integer:: comm !,nsizex,info
    comm= merge(commin,MPI_COMM_WORLD,present(commin))
    !call MPI_Comm_size( comm, nsizex, info ); write(*,*) 'mmmmmmmyyyy 1111 mpisizexxxxxx=',nsizex
    if(master_mpi) write(stdo,"(a)")'m_lmfinit: '//trim(prgnam)
    
    ! ConvertCtrl2CtrlpByPython: block
    !   use m_args,only: argall
    !   character(512):: cmdl
    !   logical:: fileexist
    !   if(master_mpi) then
    !      inquire(file='ctrl.'//trim(sname),exist=fileexist)
    !      if(.NOT.fileexist) call rx("No ctrl file found!! ctrl."//trim(sname))
    !      open(newunit=ifi,file='save.'//trim(sname),position='append')
    !      write(ifi,"(a)")'Start '//trim(prgnam)//trim(argall)
    !      close(ifi)
    !      cmdl=trim(cmdpath)//'ctrl2ctrlp.py '//trim(argall)//'<ctrl.'//trim(sname)//' >ctrlp.'//trim(sname)
    !      write(stdo,"(a)")'cmdl for python='//trim(cmdl)
    !      call system(cmdl) !See  results ctrlp.* given by ctrl2ctrl.py 
    !   endif
    !   call MPI_BARRIER( comm, ierr)
    ! end block ConvertCtrl2CtrlpByPython

    ReadCtrlp: block ! Readin ctrlp
      character(10000):: recrdx=''
      integer::ixx,lenmax
      open(newunit=ifi,file='ctrlp.'//trim(sname))       !read(ifi,*) nrecs2,reclnr !nrecs,
      lenmax=-9999
      ixx=0
      do 
         read(ifi,"(a)",end=1011)recrdx
         ixx=ixx+1
         lenmax=max(len_trim(recrdx),lenmax)
      enddo
1011  continue
      reclnr=lenmax !max record size
      nrecs2=ixx    !number of records
      rewind ifi
      allocate(character(reclnr):: recrd(nrecs2))
      do i = 1, nrecs2
         read(ifi,"(a)") recrd(i)
      enddo
      close(ifi)
    endblock ReadCtrlp
    Stage1GetCatok: block ! Readin Category-Token-Subtoken from recrd by rval2
      logical:: cmdopt0,cmdopt2,parmxp
      integer:: iprint,isw,ncp,nmix,broy,n,n1,n2,n3
      real(8):: avwsr,rydberg,wt(2),beta
      character(8):: fnam,xn
      call gtv2_setrcd(recrd)
      call rval2('STRUC_NSPEC',rr=rr, nreq=1);  nspec=nint(rr)!number of specs
      call rval2('STRUC_NBAS', rr=rr, nreq=1);  nbas=nint(rr) !number of sites  sites=>class=>specs 
      call rval2('HAM_NSPIN',  rr=rr, defa=[real(8):: 1]);  nsp=nint(rr) 
      allocate(pnuall(n0,nsp,nbas),pnzall(n0,nsp,nbas),pnusp(n0,nsp,nspec),qnu(n0,nsp,nspec), pzsp(n0,nsp,nspec),idmod(n0,nspec), &
           rsmh1(n0,nspec),eh1(n0,nspec),rsmh2(n0,nspec),eh2(n0,nspec), rg(nspec),rsma(nspec),rfoca(nspec),&
           rmt(nspec), spec_a(nspec),z(nspec),nr(nspec),eref(nspec), coreh(nspec),coreq(2,nspec), &
           idxdn(n0,nkap0,nspec), idu(4,nspec),uh(4,nspec),jh(4,nspec), &
           cstrmx(nspec),frzwfa(nspec), kmxt(nspec),lfoca(nspec),lmxl(nspec),lmxa(nspec),&
           lmxb(nspec),nmcore(nspec),rs3(nspec),eh3(nspec),&
           lpz(nspec),lpzex(nspec),nkapii(nspec),nkaphh(nspec),slabl(nspec),&
           pos(3,nbas),ispec(nbas),ifrlx(3,nbas),iantiferro(nbas))
      idu=0; uh=0d0; jh=0d0; rs3=0.5d0; eh3=0.5d0; pnusp=0d0; pzsp=0d0; qnu=0d0; lpz=0; lpzex=0; cstrmx=F 
      nkapii=1; nkapi=1; rsmh1 = 0d0; rsmh2 = 0d0; eh1  = 0d0; eh2 = 0d0; idmod=0; rfoca = 0d0; rg=0d0
      call rval2('IO_VERBOS' , rr=rr, defa=[real(8)::  30]); verbos=nint(rr)
      call rval2('IO_TIM'    , rr=rr, defa=[real(8)::  1 ]); io_tim=nint(rr)
      call rval2('STRUC_ALAT', rr=rr, nout=n);  alat=rr  !   lattice parameter, in a.u.
      call rval2('STRUC_DALAT',rr=rr, nout=n);  dalat=rr !adding to ALAT
      call rval2('STRUC_PLAT', rv=rv, nreq=9);  plat=reshape(rv,shape(plat))
      call rval2('OPTIONS_HF' ,rr=rr, defa=[real(8):: 0]);  lhf= nint(rr)==1 ! for non-self-consistent Harris'
      call rval2('HAM_REL',    rr=rr, defa=[real(8):: 1]);  lrel=nint(rr)!lrel=0 non relativistic. lrel=1 scalar relativistic
      call rval2('HAM_SO',     rr=rr, defa=[real(8):: 0]);  lso=nint(rr)
      call rval2('HAM_SOCAXIS',rv=rv, defa=[0d0,0d0,1d0]);  socaxis=rv
      ! SOnote: In general cases (except 001), we nees spin-symmetirc radial functions because we use  <up|Lz|up> for the place of <up|Lz|dn>, for example.
      ! If we like to take into account spin-dependent radial functions, we need to calculate ohsozz%soffd, ohsopm%sdiag.
      !  See 'call gaugm' at the end of augmat.
      call rval2('HAM_GMAX',  rr=rr, defa=[0d0]);          lat_gmaxin=rr
      call rval2('HAM_FTMESH',rv=rv, defa=[0d0,0d0,0d0]);  ftmesh=rv
      call rval2('HAM_TOL',   rr=rr, defa=[1d-6]); tolft=rr  
      call rval2('HAM_FRZWF', rr=rr, defa=[real(8):: 0]); ham_frzwf= nint(rr)==1 
      call rval2('HAM_XCFUN', rr=rr, defa=[real(8):: 2]); lxcf=nint(rr) ! XC functional. 1:Ceperly-Alder, 2:Barth-Hedin, 103:PBE-GGA
      call rval2('HAM_FORCES',rr=rr, defa=[real(8):: 0]); lfrce=nint(rr) !0 or 1
      call rval2('HAM_RDSIG', rr=rr, defa=[real(8):: 1]); lrsigx=nint(rr)
      call rval2('HAM_ScaledSigma', rr=rr, defa=[1d0]   ); scaledsigma=rr
      call rval2('HAM_EWALD', rr=rr, defa=[real(8):: 0]); ham_ewald= nint(rr)==1
      call rval2('HAM_OVEPS', rr=rr, defa=[1d-7]);   oveps=rr
      call rval2('HAM_PWMODE',rr=rr, defa=[real(8):: 0]);  pwmode=nint(rr)
      call rval2('HAM_PWEMAX',rr=rr, defa=[real(8):: 0]);  pwemax=rr
      call rval2('HAM_READP', rr=rr, defa=[real(8):: 0]); readpnu= nint(rr)==1
      call rval2('HAM_READPSKIPF', rr=rr, defa=[real(8):: 1]); readpnuskipf= nint(rr)==1
      call rval2('HAM_V0FIX', rr=rr, defa=[real(8):: 0]); v0fix =  nint(rr)==1
      call rval2('HAM_PNUFIX',rr=rr,defa=[real(8):: 0]); pnufix=  nint(rr)==1
      avw = avwsr(plat,alat,vol,nbas) !vol cell volume
      specloop: do j=1,nspec !SPEC_ATOM j is spec index. In SPEC category, we do j=j+1 after we find ATOM=xx. See ctrl2ctrlp.py
         if(master_mpi) write(stdo,"(a,g0)")'=== SPEC =',j
         call rval2('SPEC_ATOM@'//xn(j),ch=ch); slabl(j)=trim(adjustl(ch))
         call rval2('SPEC_Z@'   //xn(j),rr=rr); z(j)=rr
         call rval2('SPEC_R@'   //xn(j),rr=rr, nout=n1); if(n1==1) rmt(j)=rr
         call rval2('SPEC_R/W@' //xn(j),rr=rr, nout=n2); if(n2==1.and.n1==0) rmt(j)=rr*avw
         call rval2('SPEC_R/A@' //xn(j),rr=rr, nout=n3); if(n3==1.and.n2==0.and.n1==0) rmt(j)=rr*alat
         call rval2('SPEC_A@'   //xn(j),rr=rr, defa =[0d0]);        spec_a(j)=rr
         call rval2('SPEC_NR@'  //xn(j),rr=rr, defa=[real(8)::0]);  nr(j)=nint(rr)
         call rval2('SPEC_RSMH@'//xn(j),rv=rv,nout=nnx); rsmh1(1:nnx,j)=rv       ! 1st MTO sets
         nn1 = findloc([(rsmh1(i,j)>0d0,i=1,nnx)],value=.true.,dim=1,back=.true.)! 1st MTO sets
         call rval2('SPEC_EH@'  //xn(j),rv=rv,nreq=nn1); eh1(1:nn1,j)=rv         ! 1st MTO sets 
         call rval2('SPEC_RSMH2@'//xn(j),rv=rv,nout=nnx); rsmh2(1:nnx,j)=rv      ! 2nd MTO sets 
         nn2 = findloc([(rsmh2(i,j)>0d0,i=1,nnx)],value=.true.,dim=1,back=.true.)! 2nd MTO sets 
         call rval2('SPEC_EH2@' //xn(j),rv=rv, nreq=nn2); eh2(1:nn2,j)=rv        ! 2nd MTO sets 
         if(nn2>0) nkapi=2
         if(nn2>0) nkapii(j)=2
         call rval2('SPEC_LMX@'//xn(j), rr=rr, defa=[real(8):: 999]);     lmxb(j)=min(nint(rr), max(nn1,nn2)-1)! lmax for MTO basis
         call rval2('SPEC_LMXA@'//xn(j),rr=rr, defa=[real(8):: lmxb(j)]); lmxa(j)=nint(rr)
         if(rmt(j)==0d0) lmxa(j)=-1                                                        ! lmax for augmentation
         call rval2('SPEC_LMXL@'//xn(j),rr=rr, defa=[real(8):: lmxa(j)]); lmxl(j)=nint(rr) ! lmax for accumulating rho,V in sphere
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
         ! 2019 Automatic turn off lda+u mode by +10 ; Enforce uh=jh=0 when sigm exist !!
         call rval2('SPEC_UH@'//xn(j),  rv=rv,nout=n); uh(1:n,j)=rv ! Hubbard U for LDA+U
         call rval2('SPEC_JH@'//xn(j),  rv=rv,nout=n); jh(1:n,j)=rv ! Exchange parameter J for LDA+U
         call rval2('SPEC_C-HOLE@'//xn(j),ch=ch); coreh(j)=trim(adjustl(ch)) ! Channel for core hole
         call rval2('SPEC_C-HQ@'//xn(j),rv=rv,defa=[-1d0,0d0]) ; coreq(:,j)=rv
         !   NOTE: 'Charge and Moment of core hole:'// &
         !          Q(spin1) = full + C-HQ(1)/2 + C-HQ(2)/2 
         !          Q(spin2) = full + C-HQ(1)/2 - C-HQ(2)/2
         call rval2('SPEC_EREF'//xn(j), rr=rr,defa=[0d0]); eref(j)=rr
      enddo specloop
      lmxbx=maxval(lmxb)
      ibasloop: do j = 1, nbas
         if(master_mpi) write(stdo,"(a,g0)")'=== SITE =',j
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
      call rval2('BZ_NKABC',rv=rv, nout=n); bz_nabcin(1:n)=nint(rv) !No. q points along 3 lattice vectors.=no. divisions for F.T. mesh
      call fill3in(n,bz_nabcin) !filled to the end if n<3
      call rval2('BZ_BZJOB',rv=rv, nout=n); bz_lshft(1:n)=nint(rv)! '=0 centers BZ mesh at origin, =1 centers off origin' 
      call fill3in(n,bz_lshft)
      call rval2('BZ_METAL', rr=rr, defa=[real(8):: 3]); bz_lmet=nint(rr) !0 or 3. '0 insulator only; 3 for metal
      call rval2('BZ_TETRA', rr=rr, defa=[real(8):: 1]); bz_tetrahedron= nint(rr)==1 ! tetrahedron switch 1 or 0.
      if(cmdopt0('--tdos') .OR. cmdopt0('--pdos')) then
         write(stdo,*)' --tdos or --pdos enforces BZ_METAL=3 and BZ_TETRA=T'
         bz_lmet=3
         bz_tetrahedron=.true.
      endif
      call rval2('BZ_N', rr=rr, defa=[real(8):: 0]); bz_n=nint(rr) 
      !          BZ_N>0: Polynomial order for Methfessel-Paxton sampling', N=0: Conventional Gaussian sampling
      !          BZ_N<0: Broadening by Fermi-Dirac distribution'.     Use in conjunction with BZ_W
      call rval2('BZ_W', rr=rr, defa=[5d-3]); bz_w=rr ! For BZ_N>=0,Line broadening for sampling. For BZ_N<0, Temperature for Fermi distribution (Ry)
      call rval2('BZ_ZBAK',  rr=rr, defa=[0d0]); zbak=rr !Homogeneous background charge to shift Fermi energy.
      call rval2('BZ_SAVDOS',rr=rr, defa=[real(8):: 0]); ldos=nint(rr)! '0(F) or 1(T): Write dos.tot.* file (settings are NPTS and DOS)'
      call rval2('BZ_NPTS',  rr=rr, defa=[real(8):: 2001]); bz_ndos=nint(rr) !'No. DOS points (sampling integration)')
      call rval2('BZ_DOSMAX',rr=rr, defa=[40d0/rydberg()]); bz_dosmax=rr ! Maximum energy to which DOS accumulated, relative to Efermi
      call rval2('BZ_EFMAX', rr=rr, defa=[5d0]); bz_efmax=rr !Find evecs up to efmax' 
      call rval2('BZ_NEVMX', rr=rr, defa=[real(8):: 0]); bz_nevmx=nint(rr)! Find at most nevmx eigenvec. NEVMX=0:internal default, NEVMX<0:no eigenvecs
      if( cmdopt0('--tdos') .OR. cmdopt0('--pdos') .OR. cmdopt0('--zmel0')) bz_nevmx=999999
      call rval2('BZ_FSMOM',rr=rr, defa=[NULLR]); bz_fsmom=rr !Size of Fixed-spin moment (fixed-spin moment method tured on when BZ_FSMOM/=NULLR)
      call rval2('BZ_FSMOMMETHOD', rr=rr, defa=[real(8):: 0]); bz_fsmommethod=nint(rr) !'Method of Fixed-spin moment 0:original 1:discrete')
      call rval2('SYMGRP', ch=ch); symg=adjustl(ch) ! Generators for symmetry group'
      call rval2('SYMGRPAF', ch=ch); symgaf=adjustl(ch) ! Extra Generator for adding anti ferro symmetry'
      afsym=merge(.true.,.false., len_trim(symgaf)/=0.and.lso/=1)
      
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
      !            Anderson: A[nmix][,b=beta][,bv=betv][,w=w1,w2][,nam=fn][,k=nkill]
      !            Broyden:  B[nmix][,b=beta][,bv=betv][,wc=wc][,w=w1,w2][,nam=fn][,k=nkill]
      ch=trim(adjustl(ch))//' -1' ! nmix=-1 default
      broyinit=merge(1,0,ch(1:1)=='B'.or.ch(1:1)=='b') !=0 Anderson,=1 broyden
      read(ch(2:),*) nmixinit
      call rval2('ITER_CONV', rr=rr,defa=[1d-4]); etol=rr   !Tolerance in energy change from prior iteration for self-consistency')
      call rval2('ITER_CONVC',rr=rr,defa=[1d-4]); qtol=rr   !Tolerance in output-input charge for self-consistency')
      call rval2('ITER_UMIX',rr=rr,defa=[.5d0]) ; mix_umix=rr !Mixing parameter for densmat in LDA+U') !2022mar9 default umix=0.5
      call rval2('ITER_TOLU',rr=rr,defa=[0d0])  ; mix_tolu=rr !Tolerance for densmat in LDA+U'
      call rval2('ITER_b',   rr=rr,defa=[1d0])  ; betainit=rr
      if(betainit/=1d0) bexist=.true. !out
      call rval2('ITER_wc',  rr=rr,defa=[-1d0]) ; wc=rr
      call rval2('ITER_w',   rv=rv,defa=[1d0,1d0]); wtinit=rv
      call rval2('ITER_k',   rr=rr,defa=[real(8):: -1]); killj=nint(rr)
      call rval2('DYN_MODE', rr=rr,defa=[real(8):: 0]); lrlxr=nint(rr) ! Dynamics (only for relaxation now).
      !  lrlxr=0: no relaxation, 4:relaxation(conjugate gradients), 5:relaxation(Fletcher-Powell), 6:relaxation(Broyden)
      if(lrlxr/=0) lfrce=1 !force calculation turned on
      call rval2('DYN_NIT',rr=rr,  defa=[real(8):: 1]); nitrlx=nint(rr) !'maximum number of relaxation steps (statics)'//' or time steps (dynamics)')
      call rval2('DYN_HESS',rr=rr, defa=[real(8):: 1]); rdhessr= nint(rr)==1 !'Read hessian matrix')
      call rval2('DYN_XTOL',rr=rr, defa=[1d-3]);xtolr=rr !Convergence criterion in displacements XTOL>0: use length; <0: use max val; =0: do not use')
      call rval2('DYN_GTOL',rr=rr, defa=[0d0]); gtolr=rr !Convergence criterion in gradients'GTOL>0: use length;  <0: use max val;  =0: do not use')
      call rval2('DYN_STEP',rr=rr, defa=[0.015d0]); stepr=rr !Initial (and maximum) step length'
      call rval2('DYN_NKILL',rr=rr,defa=[real(8):: 0]);nkillr=nint(rr)!'Remove hessian after NKILL iter')
      if(master_mpi)write(stdo,ftox)'mixing param: A/B nmix wt=',broyinit,nmixinit,ftof(wtinit),'beta wc killj=',&
           ftof(betainit),ftof(wc),killj
    endblock Stage1GetCatok
    Stage2SetModuleParameters: block
      integer:: isw,iprint
      logical:: cmdopt0,cmdopt2
      allocate(sspec(nspec)) !NOTE: this is in module m_fatom.f90
      call setpr0(30)       !Set initial verbose for console output.
      if    (cmdopt2('--pr=',outs))then; read(outs,*) verbos
      elseif(cmdopt2('--pr',outs)) then; read(outs,*) verbos
      elseif(cmdopt2('-pr',outs))  then; read(outs,*) verbos
      endif
      call setpr0(verbos) !Set initial verbos
      if( .NOT. master_mpi) call setpr0(-100) !iprint()=0 except master
      io_tim=0
      if(cmdopt2('--time',outs)) then ! Timing: Turns CPU timing log #1:tree depth #2:CPU times as routines execute.
         outs=trim(outs(2:))//' 999' ! --time=#1,#2     
         read(outs,*)io_tim(1),io_tim(2)
         if(io_tim(1)==999) io_tim=[5,2]
      endif
      call tcinit(io_tim(2),io_tim(1),levelinit) !Start tcn (timing monitor) 
      call tcn('m_lmfinit') 
      lcd4 = merge(T,F,prgnam/='LMFA')
      if(sum(abs(socaxis-[0d0,0d0,1d0])) >1d-6 .AND. (.NOT.cmdopt0('--phispinsym'))) &
           call rx('We need --phispinsym for SO=1 and HAM_SOCAXIS/=001. Need check if you dislike --phispinsym.')
      !TK found --phispinsym (the same radialfunctions for both spins) caused a problem to determine required number of nodes for NiO(LDA). 2023sep.
      if(cmdopt0('--zmel0')) OVEPS=0d0 !for epsmode ok?
      if(pwmode==10) pwmode=0   !takao added. corrected Sep2011
      if(prgnam=='LMFGWD') pwmode=10+ mod(pwmode,10) !lmfgw mode use 
      if(iprint()>0) write(stdo,ftox) ' ===> for --jobgw, pwmode is switched to be ',pwmode
      inquire(file='sigm.'//trim(sname),exist=sexist)
      nspecloop0: do j = 1, nspec ! Radial mesh parameters: determine default value of a
         i0 = NULLI
         xxx = NULLR
         call pshpr(0)
         call rmesh(z(j),rmt(j),lrel,.false.,nrmx,xxx,i0)
         if (xxx == .03d0) xxx = .015d0 !.025d0 jun2012 .025 to .015 as default.
         if(spec_a(j)==0d0) spec_a(j)=xxx
         i0 = 0
         call rmesh(z(j),rmt(j),lrel,.false.,nrmx,spec_a(j),i0)
         call poppr
         rg(j)   = 0.25d0*rmt(j)
         rfoca(j)= 0.4d0*rmt(j)
         if(nr(j) == 0) nr(j) = i0
         call rxx( nr(j)         > nrmx,'m_lmfinit: increase nrmx')
         call rxx((lmxl(j)+1)**2 > nlmx,'m_lmfinit: increase nlmx')!nlml
         nnx = findloc([(pzsp(i,1,j)>0d0,i=1,n0)],value=.true.,back=.true.,dim=1)
         lmxb(j) = max(lmxb(j),nnx-1) ! lmxb corrected by pzsp
         if(nnx>0) lpz(j)=1
         if(nnx>0.and.sum(floor(pzsp(1:lmxa(j)+1,1,j)/10))>0 ) lpzex(j)=1
         if(maxval(pzsp(1:n0,1,j))>10d0) lpztail= .TRUE. ! PZ +10 mode exist or not.
         nkaphh(j) = nkapii(j) + lpz(j) !number of radial basis of MTOs for j.
      enddo nspecloop0
      PnuQnuSetting: do 1111 j = 1, nspec ! Radial mesh parameters: determine default value of a
         if (lmxa(j)==-1) cycle
         !Pnu is fractional quantum number. For example, P=4.56; 4 is principle quantum number (nodenum-l);
         !    .56 is for log derivative (+inf to -inf is mapped to 0 to 1) at MT to determine radial functions.
         ! P = PrincipleQnum - 0.5*atan(dphidr/phi)/pi
         ! === Reset default P,Q in absence of explicit specification ====
         !     !  qnu is set by default p. --
         !     ! This looks too complicated. Fix this in future.
         !     ! In anyway, we expect pnusp and qnu are correctly returned (qnu does not care value of given P).
         !     ! set default pnusp. See the following section 'correct qnu'
         !     ! isp=1 means charge. isp=2 means mmom
         ReadDefaultPnuQnu: block
           if(allocated(pnuspdefault)) deallocate(pnuspdefault,qnudefault,qnudummy)
           allocate(pnuspdefault(n0,nsp),qnudefault(n0,nsp),qnudummy(n0,nsp),source=0d0)
           iqnu = merge(0,1,sum(abs(qnu(:,1,j)))<1d-8) !check initial Q is given or not.
           call defpq(z(j),lmxa(j),1,pnuspdefault,qnudefault)! qnu is given here for default pnusp.
           call defpq(z(j),lmxa(j),1,pnusp(1,1,j),qnudummy)  ! Set pnusp. 
           if(iqnu==0) qnu(:,1,j)  = qnudefault(:,1) !charge (qnu(:,:,2) is maga moment.)
           if(nsp==2) pnusp(1:n0,2,j)= pnusp(1:n0,1,j)
           if(nsp==2) pzsp (1:n0,2,j)= pzsp (1:n0,1,j)
         endblock ReadDefaultPnuQnu
         ReadPnuOverwritepnu:block 
           integer:: ifipnu,lr,iz,nspr,lrmx,isp,ispx
           real(8):: pnur,pzav(n0),pnav(n0),pzsp_r(n0,nsp,nspec),pnusp_r(n0,nsp,nspec)
           character(8):: charext
           if(prgnam /= 'LMFA'.and.ReadPnu) then
              pzsp_r =0d0
              pnusp_r=0d0
              open(newunit=ifipnu,file='atmpnu.'//trim(charext(j))//'.'//trim(sname))
              if(master_mpi) write(stdo,*)'READP=T: read pnu from atmpnu.*'
              do
                 read(ifipnu,*,end=1015) pnur,iz,lr,isp
                 if(iz==1) pzsp_r (lr+1,isp,j)= pnur ! +10d0 caused probelm for 3P of Fe.
                 if(iz==0) pnusp_r(lr+1,isp,j)= pnur
                 lrmx=lr
                 nspr=isp
              enddo
1015          continue
              pzav(1:lrmx+1)=sum(pzsp_r (1:lrmx+1,1:nspr,j), dim=2)/nspr !spin averaged
              pnav(1:lrmx+1)=sum(pnusp_r(1:lrmx+1,1:nspr,j), dim=2)/nspr
              do l=1,lrmx+1
                 if(pzav(l)>pnav(l)) pzav(l)=floor(pzav(l))+.5d0 ! push up p =floor(p)+0.5 if we have lower orbital
                 if(pzav(l)>1d-8.and.pnav(l)>pzav(l)) pnav(l)=floor(pnav(l))+.5d0
              enddo
              do ispx=1,nspr
                 do lr=1,lrmx+1
                    pzsp(lr, ispx,j) = pzav(lr)
                    if(lr>3.and.readpnuskipf)cycle
                    pnusp(lr,ispx,j)=  pnav(lr)
                 enddo
              enddo
              close(ifipnu)
           endif
         endblock ReadPnuOverwritepnu
         ! We have four cases: run >lmfa |grep conf
         !     P=Pdefault      ! qnu
         !     Pdefault < P    ! Default core is treated as valence with Pdefault.
         !     Pz < P=Pdefault ! qnu + 2*(2l+1)
         !     Pz=Pdefault < P ! qnu
         if(iqnu==0) then
            do lx=0,lmxa(j)     !correct valence number of electrons.
               if(pzsp(lx+1,1,j)<1d-8) then ! PZSP(local orbital) not exist
                  if(int(pnuspdefault(lx+1,1))<int(pnusp(lx+1,1,j))) qnu(lx+1,1,j)=0d0 !pnuspdefault is filled and no q for pnusp. (core hole case or so)
               else           !PZ exist   !     print *,'qnu=',lx,qnu(lx+1,1,j)
                  if( mod(int(pzsp(lx+1,1,j)),10)<int(pnuspdefault(lx+1,1)) ) qnu(lx+1,1,j)= qnu(lx+1,1,j)+ 2d0*(2d0*lx+1d0)
               endif
            enddo
         endif
1111  enddo PnuQnuSetting
      SkipLDAUwhenSigmExistFORidu10plus: if(sexist) then
         do j=1,nspec
            if(sum(abs(idu(:,j)))/=0) then
               do lxxx=0+1,3+1
                  if(idu(lxxx,j)>10) then
                     if(master_mpi) write(stdo,"(a,2i4)")'For IDU>10 with sigm.*, we set UH=JH=0 for l,ibas=',lxxx,j
                     uh(lxxx,j) = 0d0
                     jh(lxxx,j) = 0d0
                  endif
               enddo
               idu(0+1:3+1,j) = mod(idu(0+1:3+1,j),10)
            endif
         enddo
      endif SkipLDAUwhenSigmExistFORidu10plus
      lmxax = maxval(lmxa) !Maximum L-cutoff
      maxit=iter_maxit
      nlmax=(max(lmxbx,lmxax)+1)**2
      if (dalat == NULLR) dalat=0d0
      lat_alat=alat+dalat
      lat_avw=avw
      lat_nkqmx=lat_nkdmx
      lat_platin=plat
      lat_tolft=tolft
      lmxcg=10 ! set cg coefficients for lmf --- Choose dimensions for arrays
      lmxcy=12
      lnjcg = 62200 !  for (lmxcg .le. 10);  lnjcg = 62200; lnxcg = 7400
      lnxcg = 7400  !3400
      nlm=(lmxcy+1)**2
      allocate(rv_a_ocy(nlm),rv_a_ocg(lnjcg),iv_a_ojcg(lnjcg),iv_a_oidxcg(lnxcg))
      call sylmnc(rv_a_ocy , lmxcy ) ! Clebsh-Gordon coefficients for lmf part
      call scg(lmxcg , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg ) !set CG coefficients for lmf part.
      if (master_mpi) then
         inquire(file='sigm.'//trim(sname),exist=sexist)
         if (lrsigx/=0 .AND. ( .NOT. sexist) ) then
            write(stdo,*)' bndfp (warning): no sigm file found ... LDA calculation only'
            lrsigx = 0
         endif
         ham_lsig=lrsigx
      endif
      call mpibc1_int(ham_lsig,1,'bndfp_ham_lsig')
      ham_scaledsigma=scaledsigma
      ham_pwmode=pwmode
      !r  pwmode Controls PW part of basis.  
      !         0 => no PW part of basis
      !         1 => include PWs in basis by |q|^2   <pwemax
      !         2 => include only PWs in basis   
      !      11,12=>  we do q cutoff by      |q+G|^2 <pwemax
      ham_oveps=oveps
      idxdn=0 !=1  when MTO of (lp1,ik,j) included in the Hamiltonian
      do j=1,nspec
         do lp1 = 1, lmxb(j)+1
            do  ik = 1, nkapii(j)
               if(ik==1.AND. rsmh1(lp1,j)>0) idxdn(lp1,ik,j)=1
               if(ik==2.AND. rsmh2(lp1,j)>0) idxdn(lp1,ik,j)=1
            enddo
            if(pzsp(lp1,1,j)>0)              idxdn(lp1,nkaphh(j),j)=1 
         enddo
      enddo
      allocate(rsmv,source=.5d0*rmt)
      if(master_mpi) write(stdo,*)
      if(master_mpi) write(stdo,ftox)'pnu === pnu setting ==='
      if(master_mpi) write(stdo,ftox)'pnu   ibas isp   pnu(0:lmxa)  pz(0:lmxa)'
      do j=1,nbas
         is=ispec(j) 
         pnuall(:,1:nsp,j) = pnusp(1:n0,1:nsp,is)
         pnzall(:,1:nsp,j) = pzsp(1:n0,1:nsp,is)
         if(master_mpi) write(stdo,ftox)('pnu: ',j,isp,ftof(pnuall(1:lmxa(is)+1,isp,j),3),' ',&
              ftof(pnzall(1:lmxa(is)+1,isp,j),3), isp=1,nsp)
      enddo
      !! ... Suppress symmetry operations for special circumstances
      ! addinv=T only when we expect Hamiltonian is real even with keeping spin direction.
      ! addinv=T means psi* is eigenfunction.(Time-Reversal with keeping spin).
      ! 2022-jan-20 new setting of addinv (addinv =.not.ctrl_noinv)
      !Add inversion to get !When we have TR, psi_-k(r) = (psi_k(r))^* (when we have SO/=1).
      !                      density |psi_-k(r)|^2 = |psi_k^*(r)|^2

      ! # NOINV={noinv} (If F, we enfoce |psi_\sigma^fk|^2=|psi_\sigma^{-fk}|^2)
      ! #    usually judged automatically. So this is for debugging).
      ! # If NOINV=F, we assume H=H^* for each spin (k <-> -k symmetry). 
      ! #   This is the case in DFT without spin-orbit coupling.
      ! #   Then we use |psi_\sigma^fk|^2=|psi_\sigma^{-fk}|^2 to reduce comutational time.
      ! # If so/=0, lmf automatically sets NOINV=T.
      ! # You must get correct results evenif NOINV=T in the case of DFT without SO (for debug).
      ! #
      ! # WARN: Even if sigm exist, we use NOINV=F; exactly speaking, this is not correct.
      ! #       sigm without SO can cause orbital moments, that is, breaking time-reversal symmetry.
      ! # NOTE: because of inversion in space-group symmetry, we may have 
      ! #       |phi_sigm^fk|^2 = |phi_sigm^{-fk}|^2. This is not for NOINV.
    
      if( (lrlxr>=1.AND.lrlxr<=3) .OR. cmdopt0('--cls') .OR. cmdopt0('--nosym') .OR. cmdopt0('--pdos')) then
         symg = 'e'
         addinv = .false. 
!   elseif( lso==0 .and. sum(abs(idu))/=0 .and. (.not.sexist) ) then ! Add inversion means Hamiltonian is real (time-reversal).
      elseif( lso==0 .and. sum(abs(idu))==0 .and. (.not.sexist) ) then ! Add inversion means Hamiltonian is real (time-reversal). !bugfix by obata 2023-12-12
         addinv=.true. 
      else
         addinv=.false. 
      endif
      sstrnsymg=trim(symg)
      nspc = merge(2,1,lso==1) ! nspc=2 for lso=1. nspx=nsp/nspc. Haittonian is ham(1:ndham*nspc, 1:ndham*nspc, 1:nsp/nspc)
      GETorbitalindex: block ! probably too confusing.
        integer:: ib,l,lmr,ia, nnrlx,lmri,ik,nnrl,nnrli,li, iorb,is,k,iorbmto,iorbe,jorb,io,nnnx
        logical:: agree
        !o   norb  :number of orbital types for ib; see Remarks
        !o   ltab  :table of l-quantum numbers for each type
        !o   ktab  :table of energy index for each type
        !o   offl  :offl(norb) offset in h to this block of orbitals
        !o   ndim  :dimension of hamiltonian for this site
        ndimx=0
        do j=1,nspec
           nnnx=0
           do li=0,n0-1
              nnnx= nnnx + sum((2*li+1)*idxdn(li+1,:,j))
           enddo
           if(nnnx>ndimx) ndimx=nnnx
        enddo
!        ndimx=maxval([ (sum([(sum((2*li+1)*idxdn(li+1,:,j)), li=0,n0-1)]), j=1,nspec) ])
        allocate(ltabx(n00,nbas),ktabx(n00,nbas),offlx(n00,nbas),ndimxx(nbas),norbx(nbas))
        norbx=0
        ndimxx=0
        nlmto=0
        do ib=1,nbas
           is = ispec(ib)
           do  ik = 1, nkap0
              do  l = 0, n0-1
                 offlx(norbx(ib)+1,ib) = -1
                 if(idxdn(l+1,ik,is)==1) then
                    offlx(norbx(ib)+1,ib) = nlmto 
                    norbx(ib)  = norbx(ib) +1
                    ndimxx(ib) = ndimxx(ib) + 2*l+1
                    ltabx(norbx(ib),ib) = l
                    ktabx(norbx(ib),ib) = ik
                    nlmto=nlmto+2*l+1
                 endif
              enddo
           enddo
           if (norbx(ib) > n00) call rx('orbl: norb> n00')
        enddo
        norbxx=maxval(norbx)   
        norbmto=sum(norbx)
        allocate(blksx(norbxx,nbas),ntabx(norbxx,nbas)) ! requires l be consecutive and kappa index be constant in contiguous block
        !ntab  :table of upper ranges for each orbital block, !  blks sum of 2*l+1
        GETntabxblksx: block ! Marks blocks of contiguous l for which rsm and e are unchanged
          do ib=1,nbas
             iorb=1
             do while(iorb<=norbx(ib))
                do jorb=iorb,norbx(ib)
                   agree = ltabx(jorb,ib)==ltabx(iorb,ib)+(jorb-iorb) .and. ktabx(jorb,ib)==ktabx(iorb,ib)
                   if(agree) iorbe=jorb
                   if(.not.agree) exit
                enddo
                ntabx(iorb,ib)= iorbe
                ntabx(iorb+1:iorbe,ib)=0
                iorb=iorbe+1
             enddo
             blksx(1:norbx(ib),ib) =[(sum( [ (2*ltabx(io,ib)+1, io=iorb,ntabx(iorb,ib)) ] ),iorb=1,norbx(ib))]
          enddo
        endblock GETntabxblksx
        !!--- make index table :norbmto is the total number of different type of MTOs !      allocate( ibasindex(ndimham))
        allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
        iorbmto=0 !      ndimham = 0 !dimension of mto part of hamiltonian!      allocate(offH(nbas+1)) !offH looks!      offH=0
        do  ib = 1, nbas
           do  iorb = 1, norbx(ib) !(ib,irob) specify a block of MTO part Hamiltonian
              iorbmto=iorbmto+1
              ibastab(iorbmto)= ib
              ltab(iorbmto)   = ltabx(iorb,ib) !angular momentum l of (ib,iorb) block
              ktab(iorbmto)   = ktabx(iorb,ib) !radial index of (ib,iorb) block
              offl(iorbmto)   = offlx(iorb,ib) !offset to (ib,iorb) block
           enddo
        enddo
        kxx=maxval([(ktabx(1:norbx(ib),ib),ib=1,nbas)])
        lxx=maxval([(ltabx(1:norbx(ib),ib),ib=1,nbas)])
        allocate(offlrev(nbas,0:lxx,kxx))
        allocate(ib_table(nlmto),l_table(nlmto),k_table(nlmto))
        do iorb=1,norbmto 
           ibas = ibastab(iorb)
           offlrev(ibas,ltab(iorb),ktab(iorb))= offl(iorb) ! ... reverse maping of offset-index for hamiltonian
           ib_table(offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ibas
           l_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ltab(iorb)
           k_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ktab(iorb)
        enddo
      endblock GETorbitalindex
      nspx  = nsp/nspc
      nvi= sum([( (lmxa(ispec(ib))+1)**2,ib=1,nbas )])
      nvl= sum([( (lmxl(ispec(ib))+1)**2,ib=1,nbas )])
      nlmxlx = maxval([( (lmxl(ispec(ib))+1)**2,ib=1,nbas )])
      pot_nlma=nvi
      pot_nlml=nvl
      allocate(jnlml(nbas), source=[1,(1+sum( (lmxl(ispec(1:i))+1)**2 ),i=1,nbas-1)]) !offset of lmxl
      deallocate(idxdn) 
      allocate(iv_a_oips(nbas),source=[(ispec(ib), ib=1,nbas)])
      seref= sum([(eref(ispec(ib)),ib=1,nbas)])
      ham_seref= seref
      call MPI_BARRIER(comm, ierr )
      if( cmdopt0('--quit=show') ) call rx0(trim(prgnam)//' --quit=show')
    endblock Stage2SetModuleParameters
    Stage3InitialSetting :block 
      logical:: cmdopt0
      integer:: iprint
      if (cmdopt0('--etot')) then
         lfrce=0 !force calculation or not
         maxit=1 !max number of iteration for electronic structure part
      endif
      if(lhf) maxit= 1
      cc=merge(274.074d0,1d10,lrel/=0) !speed of light. scalar rel 1d10 for lrel=0
      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') !compute eigenvalues for fullmesh q points 
      lrout = 1 ! =1:evaluate output density and/or KS energy 
      leks = 1
      if (bz_nevmx == -1) then !bz_nevmx ! Maximum number of eigenvalues
         lrout = 0
         leks = 0
      endif
      if(cmdopt0('--band') .OR. fullmesh) then !Switch to plot bands at specified qp
         lfrce = 0
         plbnd = 1 !band plot mode
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
      LDApU: block 
        nlibu = 0
        lmaxu = 0
        allocate(lldau(nbas))
        lldau = 0
        if(master_mpi) write(stdo,*)
        do  ib = 1, nbas
           is = ispec(ib)
           do  lx = 0, min(lmxa(is),3)
              if (idu(lx+1,is) /= 0) then
                 if (lldau(ib) ==0) lldau(ib) = nlibu+1
                 nlibu = nlibu+1
                 lmaxu = max(lmaxu,lx)
              endif
           enddo
        enddo     
        ham_nlibu=nlibu
      endblock LDApU
      ! lhh, nkapii, nkaphh (nkaphh = nkapii(1 or 2) +1) if extented local orbital exist)
      allocate(lhh(nkap0,nspec))
      lhh=-1
      do i=1,nspec
         lmxbj = lmxb(i);  call getiout(rsmh1(1,i), lmxbj+1,lhh(1,i))
         if(nkapii(i)==2)  call getiout(rsmh2(1,i), lmxbj+1,lhh(2,i))
         if(lpz(i)==1 )    call getiout(pzsp(1,1,i),lmxbj+1,lhh(nkaphh(i),i))!lh for lo
      enddo
      ShowMTOsetting:if(master_mpi) then
         write(stdo,"('mto === MTO setting ===')")
         do i=1,nspec
            lmxbj = lmxb(i)
            write(stdo,"('mto ispec lmxb lpz nkapii nkaphh=',10i5)")i,lmxb(i),lpz(i),nkapii(i),nkaphh(i)
            write(stdo,"('mto rsmh1 ',i4,100f6.2)")i, rsmh1(1:lhh(1,i)+1,i)
            write(stdo,"('mto   eh1 ',i4,100f6.2)")i,   eh1(1:lhh(1,i)+1,i)
            if(nkapii(i)==2) write(stdo,"('mto rsmh2 ',i4,100f6.2)")i, rsmh2(1:lhh(2,i)+1,i)
            if(nkapii(i)==2) write(stdo,"('mto  eh2  ',i4,100f6.2)")i,   eh2(1:lhh(2,i)+1,i)
            if(lpz(i)==1 ) write(stdo,"('mto pz    ',i4,100f6.2)")i,    pzsp(1:lhh(nkaphh(i),i)+1,1,i)
            write(stdo,"('mto lh    ',i4,100i3)")  lhh(1:nkaphh(i),i)
         enddo
      endif ShowMTOsetting
      ForceDYNsetting: block ! Atomic position Relaxation setup (DYN mode)
        integer:: i,j,k,iprint !,ifrlx(3)
        logical:: force  ! nitrlx = num of iteration cycle for atomic relaxiation (outer loop)
        if(lrlxr==0) nitrlx=0 ! =0 --> no relaxiation. Only sc calculation for given atomic position.
        if(lrlxr>0 .AND. lrlxr<4 ) call rx('lmf not set up for MD yet')
        if(nitrlx>0) then !nitrlx >0 is for atomic position relaxiation lmf.f90
           allocate(indrx_iv(2,3*nbas))
           force = lfrce>0
           if(.NOT.force .OR. lrlxr== 0) goto 9299
           if(force) then
              j = 0
              do  i = 1, nbas
                 do  k = 1, 3
                    if(ifrlx(k,i) == 1) then
                       j = j + 1
                       indrx_iv(1:2,j) = [k,i]
                    endif
                 enddo
              enddo
           endif
           natrlx = j
           if (iprint() >= 30.and.natrlx>0) then
              if    (lrlxr == 4) then; write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (conjugate gradients) ..'
              elseif(lrlxr == 5) then; write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (Fletcher-Powell) ..'
              else;                    write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (Broyden) ..' ;   endif
              write(stdo,ftox)'   relaxing',natrlx,'variables,',nitrlx,'iterations'
              write(stdo,ftox)'   x-tol g-tol step=', xtolr,gtolr,stepr
           endif
        endif
9299    continue
      endblock ForceDYNsetting
    endblock Stage3InitialSetting
    call MPI_BARRIER(comm, ierr)
    call tcx('m_lmfinit')
  end subroutine m_lmfinit_init
  pure subroutine getiout(a,iin,iout) !a(1:iout) can be nonzero.
    intent(in):: a,iin
    intent(out):: iout
    integer:: iin,iout,i
    real(8):: a(iin)
    iout= findloc([(a(i)>0,i=1,iin)],dim=1,value=.true.,back=.true.)-1
  end subroutine getiout
  pure subroutine fill3in(nin,res)!- fills res(2) or res(2:3) if res(1) or res(1:2) are given
    intent(in):: nin
    intent(out):: res
    integer :: nin,res(3)
    if (nin==2) res(3) = res(2)
    if (nin==1) res(2:3) = res(1)
  end subroutine fill3in
  pure integer function icgi(ilm1,ilm2) !starting point of cg coefficients for ilm1,ilm2
    implicit none
    intent(in)::             ilm1,ilm2
    integer:: ilm1,ilm2,indx
    indx = (max0(ilm1,ilm2)*(max0(ilm1,ilm2)-1))/2 + min0(ilm1,ilm2)
    icgi = iv_a_oidxcg(indx)
  endfunction icgi
  pure integer function icge(ilm1,ilm2) !end point of cg coefficients for ilm1,ilm2
    implicit none
    intent(in)::             ilm1,ilm2
    integer:: ilm1,ilm2,indx
    indx = (max0(ilm1,ilm2)*(max0(ilm1,ilm2)-1))/2 + min0(ilm1,ilm2)
    icge = iv_a_oidxcg(indx+1)-1
  endfunction icge
end module m_lmfinit
