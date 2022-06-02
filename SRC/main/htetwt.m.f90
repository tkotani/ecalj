program htetwt
  !!  Calculate W-V for QSGW mode.
  !! We calculate chi0 by the follwoing three steps.
  !!  tetwt5: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part

  !      use m_readeps,only: read_eps, epsinv, w_mu, llmat2=>llmat,deallocate_eps
  !!

  use m_ReadEfermi,only: readefermi,ef
  use m_readqg,only: readngmx,readqg
  use m_hamindex,only:   Readhamindex
  use m_readeigen,only: init_readeigen,init_readeigen2,readeval
  use m_read_bzdata,only: read_bzdata, ! & <--- 'call read_bzdata' sets up following data.
  ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,qbas,ginv, &
       dq_,qbz,wbz,qibz,wibz, &
       ntetf,idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
  !     &     idteti, nstar,irk,nstbz
  use m_genallcf_v3,only: genallcf_v3, &
       nclass,natom,nspin,nl,nn,ngrp, &
       nlmto,nlnmx, nctot,niw, ! & nw_input=>nw,
  alat, delta,deltaw,esmr,symgrp,clabl,iclass, ! & diw,dw,
  invg, il, in, im, nlnm, &
       plat, pos, ecore, symgg

  use m_keyvalue,only: getkeyvalue
  use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
  use m_readqgcou,only: readqgcou

  !! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
  use m_rdpp,only: rdpp,    ! & "call rdpp" generate following data.
  nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr
  !! Generate matrix element for "call get_zmelt".
  use m_zmel,only:    ! & these data set are stored in this module, and used when
  nband,itq,ngcmx,ngpmx,    ppovlz, &
       ppbir,shtvg, miat,tiat , ntq
  !! frequency
  use m_freq,only: getfreq, &
       frhis,freq_r,freq_i, nwhis,nw_i,nw,npm !output of getfreq
  !! antiferro
  !      use m_anf,only: anfcond,
  !     & laf,ibasf !,ldima,pos,natom
  !! tetwt
  use m_tetwt,only: tetdeallocate,gettetwt, ! & followings are output of 'L871:call gettetwt')
  whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb
  !! w0 and w0i (head part at Gamma point)
  use m_w0w0i,only: w0w0i, &
       w0,w0i,llmat

  !! MPI
  use m_mpi,only: MPI__hx0fp0_rankdivider2Q,MPI__hx0fp0_rankdivider2S, &
       MPI__Qtask,MPI__InitializeQSPBM,MPI__Finalize,MPI__root, &
       MPI__Broadcast,MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,MPI__rank,MPI__size, &
       MPI__Qranktab,MPI__consoleout,MPI__Ss,MPI__Se, MPI__allreducesumS, &
       MPI__barrier, MPI__rankQ,MPI__rootQ,MPI__rootS
  !! q0p
  use m_readq0p,only: readq0p, &
       wqt,q0i,nq0i ,nq0iadd,ixyz
  use m_lgunit,only: m_lgunit_init

  implicit none
  integer,allocatable:: nwgt(:,:)
  integer::iopen,maxocc2,iclose, ixc,iqxini,iqxend, &
       ifhbe,  nprecb,mrecb,mrece,nlmtot,nqbzt,! & nband,
  i,nq0ix,ngrpmx,mxx,nqbze,nqibze,ini,ix,ngrpx ! & ngcmx,
  ,nblochpmx,ndummy1,ndummy2,ifcphi,is,nwp, ! & ifvcfpout,,mdimx,nbloch
  ifepscond,nxx,ifvxcpout,ifgb0vec &
       ,nw0,iw,ifinin,iw0,noccxv,noccx &
       ,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas &
       ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0 &
       ,kx,isf,kqxx,kp,job,nwmax ! & ,ifev1,ifev2 !,nhwtot
  ,ihis,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode &
       ,  ifchipmlog ,   nw_w,nwmin  ! ,ngpmx
  real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm, dwry,dwh,omg2, q(3),  qgbin(3),qx(3)
  real(8):: ua=1d0          ! this is a dummy.
  integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2), ndble=8
  integer,allocatable :: ngveccB(:,:), iqib(:),ifppb(:) !,lx(:) ngvecc(:,:),
  complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:), &
       zw(:,:),zw0(:,:), zxq(:,:,:),zxqi(:,:,:)
  real(8),allocatable :: eqt(:), ! & ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:)
  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:), &
       ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:), qbze(:,:),qibze(:,:) !,ecore(:,:)
  !     &  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
  complex(8),allocatable :: rcxq(:,:,:,:)
  complex(8) :: fff,img=(0d0,1d0)
  complex(8),allocatable :: wwk(:,:,:)
  real(8) ::qbzx(3)
  logical :: debug
  !      integer,allocatable:: ibasf(:)
  logical :: realomega, imagomega
  complex(8),allocatable:: zzr(:,:),ppovl(:,:),ppovlzinv(:,:) !,ppovlz(:,:)
  complex(8) :: epxxx,vcmean
  character(9) :: fileps
  character(15) :: filepsnolfc
  !      logical :: paralellx0=.true. !, hist
  character(5) :: charnum5
  character(20):: xxt
  real(8) :: Emin, Emax      ,emax2,emin2
  !     integer :: iSigma_en  !sf..21May02  !iSigma_en is integer
  ! arameter stored in GWIN_V2
  ! hich determines approximation for  self-energy.
  ! elf-energy should be made hermitian for energies to be real
  ! xx  !iSigma_en==0 SE_nn'(ef)+img integral:delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
  ! xx  !iSigma_en==1 SE_nn'(ef)+delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
  ! Sigma_en==2 [SE_nn'((e_n+e_n')/2)+h.c.]/2
  ! Sigma_en==3 [(SE_nn'(e_n)+SE_nn'(e_n'))/2+h.c.]/2
  real(8) :: omg2max,omg1max,wemax
  logical::imagonly=.false. , noq0p !,readgwinput
  integer::nwin, incwfin, verbose,nbcut,nbcut2,ifpomat,nnmx,ikpo,nn_,noo,iqxxx,nomx
  !      real(8)::efin
  logical :: nolfco=.false.
  integer:: isp1,isp2, ngc,mrecg ! bzcase,
  real(8)::  quu(3),deltaq(3),qqq(3)=0d0 !
  complex(8),allocatable:: wgt(:,:,:)
  real(8),allocatable:: qbz2(:,:)
  logical :: qbzreg
  !    logical ::smbasis !smbasis will be implemented in m_zmel.f which generates <phi|phi M>
  real(8):: q_r(3)
  complex(8),allocatable:: pomat(:,:)
  logical   :: timereversal,onceww
  integer :: jpm,ncc
  real(8) :: frr !, sciss
  integer :: ngb0,ifvcoud,idummy,igb1,igb2,ngb_in,nmbas1,nmbas2,iq0,ifisk,iqx,ig,nmbas1x !ifepstinv,
  complex(8),allocatable:: zcousq(:,:),epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:),eemat(:,:),zcousq0(:,:)
  real(8),allocatable:: vcousq(:),vcousq0(:),vcoudummy(:)
  real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq
  !! Eq.(40) in PRB81 125102
  !     complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
  !     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
  complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:), &
       w_k(:),w_ks(:),w_kI(:), w_ksI(:)
  complex(8),allocatable:: llw(:,:), llwI(:,:),aaamat(:,:)
  integer:: lxklm,nlxklm,ifidmlx,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,niwxxx,iwx,icc1,icc2
  complex(8):: vc1vc2
  integer,allocatable:: neibz(:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
  real(8),allocatable:: aik(:,:,:,:)
  integer,allocatable:: aiktimer(:,:)
  integer:: l2nl, nmbas_in , iqxendx,imb2 !iqqv,
  logical:: eibz4x0,tiii,iprintx,chipm=.false.,iqinit,localfieldcorrectionllw
  real(8)::qvv(3),ecut,ecuts,hartree,rydberg,pi
  character(128):: vcoudfile
  integer :: iqeibz
  complex(8):: epslfc, axxx(10)
  integer:: src,dest
  !      integer:: ifw0w0i
  logical :: symmetrize,eibzmode
  real(8):: schi=-9999 !dummy
  integer:: i_reduction_npm, i_reduction_nwhis,  i_reduction_nmbas2
  logical:: crpa
  integer,allocatable :: iclasst(:), invgx(:)
  integer:: ibasx,ificlass,ifile_handle,ifiq0p
  complex(8),allocatable:: ppovl_(:,:)
  logical:: tetra !,readw0w0itest=.false.
  integer::nw_ixx,nwxx

  logical:: w4pmode
  complex(8),allocatable:: wmu(:,:),wmuk(:,:)
  integer:: ifw4p,ngbq0,igb
  real(8):: qv(3,3)

  real(8)::ebmx
  integer:: nbmx,mtet(3)
  real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)

  logical:: eginit=.true.
  real(8),allocatable,save:: gfmat(:,:)
  complex(8),allocatable:: rcxqin(:)
  real(8):: egauss
  integer:: imbas1,imbas2,ipm

  !-------------------------------------------------------------------------
  call MPI__InitializeQSPBM()
  call m_lgunit_init()
  ! IME0_1001 ProgAll
  call MPI__consoleout('htetwt')
  call cputid (0)
  hartree= 2d0*rydberg()
  pi     = 4d0*datan(1d0)
  fourpi = 4d0*pi
  sqfourpi=sqrt(fourpi)
  write(6,*) ' --- htetwt ----------------'
  !! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
  ! IME0_11001 readbzdata
  call read_BZDATA()
  ! IME1_11001 "readbzdata"
  ! IME0_12001 Q0P
  !! Use regular mesh even for bzcase==2 and qbzreg()=T
  !!     off-regular mesh for bzcase==1 and qbzreg()=F
  !      if( ( bzcase()==2.and.qbzreg() )       .or.
  !     &     ( bzcase()==1.and.(.not.qbzreg()))      ) then
  !! this mechanism for qbzreg=F is too complicated. We may need to modify difinition of qbz for qbzreg=F.
  if( .NOT. qbzreg()) then ! set off-gamma mesh
     deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
     do i=1,nqbz
        qbz(:,i) = qbz(:,i) - deltaq/2d0
        write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
     enddo
  endif
  write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp

  !! === Readin by genallcf ===
  !! See "use m_genallcf_v3" at the begining of this routine
  !! We set basic data.
  incwfin= 0                !use ForX0 for core in GWIN
  !--- EFERMI
  call readefermi()
  call genallcf_v3(incwfin) !in module m_genallcf_v3
  tpioa=2d0*pi/alat
  if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
  debug=.false.
  if(verbose()>=100) debug= .TRUE. 
  if(debug) write(6,*)' end of genallc'
!!!!  WE ASSUME iclass(iatom)= iatom,  nclass = natom.  !!!!!!!!!!!!!!!!!!!!!!!!!
  if(nclass /= natom) call rx( ' htetwt: nclass /= natom ')

  !! --- Readin Offset Gamma --------
  call readq0p()
  write(6,"(' ### nqibz nq0i nq0iadd=', 3i5)")nqibz,nq0i,nq0iadd
  nqbze  = nqbz *(1 + nq0i+nq0iadd)
  nqibze = nqibz + nq0i+nq0iadd
  allocate( qbze(3, nqbze), qibze(3, nqibze))
  qbze(:,1:nqbz)  = qbz(:,1:nqbz)     ! call dcopy(3*nqbz, qbz,  1, qbze,1)
  qibze(:,1:nqibz)= qibz(:,1:nqibz) ! call dcopy(3*nqibz,qibz, 1, qibze,1)
  do i = 1,nq0i+nq0iadd
     qibze(:,nqibz+i)  = q0i(:,i)
     ini = nqbz*(1 + i -1)
     do ix=1,nqbz
        qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
     enddo
  enddo
  iqxend = nqibz + nq0i + nq0iadd
  write(6,*) ' nqibz nqibze=',nqibz,nqibze
  ! IME1_13001 "mptauof"
  ! IME0_14001 init_readeigen
  call Readhamindex()
  call init_readeigen()

  !! We get frhis,freq_r,freq_i, nwhis,nw,npm  by getfreq
  realomega = .true.
  imagomega = .true.
  tetra     = .true.
  call findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
  if( .NOT. qbzreg()) then
     allocate(qbz2(3,nqbz))
     do iq=1,nqbz
        qbz2(:,iq)=qbz(:,iq)+dq_
     enddo
     call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
     emax=max(emax,emax2)
     emin=min(emin,emin2)
     deallocate(qbz2)
  endif
  if (nctot > 0) Emin=minval(ecore(:,1:nspin))
  omg2max = (Emax-Emin)*.5d0+.2d0  ! (in Hartree) covers all relevant omega, +.2 for margin
  if(MPI__root) write(6,"(' emin emax omega2max=',3f13.5)") emin, emax, omg2max
  call getwemax(.true.,wemax) !wemax is to determine nw !real axis divisions
  if(MPI__root) write(6,"(' wemax=  ',f13.4)") wemax
  call getfreq(.false.,realomega,imagomega,omg2max,wemax,niw,ua,MPI__root)!tetra,
  !! Write freq_r
  if(realomega .AND. mpi__root) then
     open(newunit=ifif,file='freq_r') !write number of frequency points nwp and frequensies in 'freq_r' file
     write(ifif,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
     do iw= nw_i,-1
        write(ifif,"(d23.15,2x,i6)") -freq_r(-iw),iw
     enddo
     do iw= 0,nw
        write(ifif,"(d23.15,2x,i6)") freq_r(iw),iw
     enddo
     close(ifif)
  endif
  !!
  nspinmx = nspin
  nwp = nw+1
  iqxini=1
  eibzmode = eibz4x0()

!!! nov2016 moved from tetwt5 --> here
  !      call getkeyvalue("GWinput","nband_chi0",nbmx, default=nband )
  !      call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
  mtet=(/1,1,1/)
  !      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
  ! multitet=T ==> micro tetrahedron method (divided-tetrahedron). Not used so much now...
  allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz))

  !! === Use of symmetry. EIBZ procedure PRB81,125102 ===
  !!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
  if(eibzmode) then
     !! commentout block inversion Use iqxendx=iqxend because of full inversion
     iqxendx=iqxend
     allocate( nwgt(nqbz,iqxini:iqxendx), ! & qeibz(3,nqbz,iqxini:nqibz),neibz(iqxini:nqibz),
     igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx), &
          eibzsym(ngrp,-1:1,iqxini:iqxendx))
     iprintx=.false.
     write(6,*)
     write(6,"('=== Goto eibzgen === TimeRevesal switch =',l1)")timereversal()
     if(MPI__root) iprintx= .TRUE. 
     call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz, &
          timereversal(),ginv,iprintx, &
          nwgt,igx,igxt,eibzsym,tiii)
     write(6,"('Used timeRevesal for EIBZ = ',l1)") tiii
     call cputid(0)
     !!  call Spacegrouprot(symgg,ngrp,plat,natom,pos) ! all inputs.
  else !dummy allocation to overlaid -check bound !sep2014
     iqxendx=iqxend
     allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx) &
          ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
  endif
  !! rank divider
  call MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
  call MPI__hx0fp0_rankdivider2S(nspinmx)

  !! for w4phonon. all nodes have wmu array.
  w4pmode=.false.
  if(sum(ixyz)/=0) w4pmode= .TRUE. 

  !! == Calculate tetrahedron weight for x0(q,iw) == main loop 1001 for iq.
  iqinit=.true.
  write(6,'("irank=",i5," allocated(MPI__qtask)=",L5)')MPI__rank,allocated(MPI__qtask)
  do iq = iqxini,iqxend
     if(MPI__qtask(iq)) write(6,'("irank iq=",i5,i5)') MPI__rank,iq
  enddo
  !! Get ngbq0 (for q=0) and broadcast for w4p
  if( MPI__root .AND. w4pmode ) then
     q = (/0d0,0d0,0d0/)
  endif
  ! IME0_170001 do1001
  do 1001 iq = iqxini,iqxend
     if( .NOT. MPI__Qtask(iq) ) cycle
     call cputid (0)
     q = qibze(:,iq)
     if( iq==1 ) then        ! *sanity check
        if(sum(q**2)>1d-10) call rx( ' htetwt: sanity check. |q(iqx)| /= 0')
     endif
     !! ---------------------------------------------------------------
     !! === loop over spin=== =========================================
     !! ---------------------------------------------------------------
     ! IME0_180001 Do1003
     do 1003 is = MPI__Ss,MPI__Se
        write(6,"(' ### ',2i4,' out of nqibz+n0qi+nq0iadd nsp=',2i4,' ### ')") &
             iq, is, nqibz + nq0i+nq0iadd, nspin
        if(debug) write(6,*)' niw nw=',niw,nw
        isf = is
        !! Tetrahedron weight.
        !! output
        !!     nbnbx
        !!     ihw(ibjb,kx): omega index, to specify the section of the histogram.
        !!     nhw(ibjb,kx): the number of histogram sections
        !!     jhw(ibjb,kx): pointer to whw
        !!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
        !!     : histogram weights for given ib,jb,kx for histogram sections
        !!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
        !            write(6,*) ' --- goto x0kf_v4hz ---- newaniso= ',newaniso2
        !! input
        !!     ekxx1 for   rk,is
        !!     ekxx2 for q+rk,isf
        do kx = 1, nqbz
           call readeval(qbz(:,kx),   is,  ekxx1(1:nband, kx) )
           call readeval(q+qbz(:,kx), isf, ekxx2(1:nband, kx) )
        enddo
        call gettetwt(q,iq,is,isf,nwgt(:,iq),ekxx1,ekxx2,eibzmode)
        !$$$          call gettetwt(q,iq,is,isf,nwgt(:,iq),frhis,nwhis,npm,
        !$$$     i     qbas,ginv, ef, nqibz, nband,ekxx1,ekxx2, nctot,ecore,
        !$$$     i     nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
        !$$$     i     nbmx,ebmx,mtet,eibzmode) !nov2016
1003 enddo
     write(6,*) 'end of spin-loop nwp=',nwp !end of spin-loop
     ! IME1_180001 Do1003
1001 enddo
  call cputid(0)
  write(6,*) '--- end of htetwt --- irank=',MPI__rank
  call flush(6)
  !      call MPI__Finalize
  call rx0('OK! --- end of htetwt ---')
end program htetwt
