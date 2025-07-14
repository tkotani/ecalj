module m_hvccfp0
  use m_ftox
  use m_lgunit,only:stdo
  use m_cputm,only:cputm
  public hvccfp0
  private
  contains
subroutine hvccfp0() bind(C)  ! Coulomb matrix. <f_i | v| f_j>_q.  ! output  VCCFP : the coulomb matrix vcoul(nblochpmx,nblochpmx) for all qibz.
  !    strx: structure constant for e=0 (means 1/|r-r'| )
  use m_xlgen,only:lgen
  use m_lattic,only:lctoff
  use m_genallcf_v3,only: Genallcf_v3,  alat,nbas=>natom, bas=>pos
  use m_hamindex,only:    Readhamindex,plat,qlat
  use m_hamindex0,only:    Readhamindex0
  use m_read_bzdata,only: Read_bzdata, ginv,nqbz,qbz,nqibz,qibz,nq0i,wqt=>wt,q0i,nq0iadd
  use m_readqg,only:   readqg,readngmx
  use m_mpi,only: MPI__Initialize,mpi__root, MPI__Broadcast,mpi__rank,mpi__size,MPI__consoleout,ipr,comm !,mpi__iend,mpi__iini !,mpi__getrange
  use m_readgwinput,only: ReadGwinputKeys, keeppositivecou
  use m_lgunit,only: m_lgunit_init
  use m_vcoulq,only: vcoulq_4,mkjb_4,mkjp_4,genjh
  use m_pwmat,only: mkppovl2
  use m_hvccfp0_util,only: mkb0,strxq
  use m_nvfortran,only:findloc
  use m_gpu,only: gpu_init
  implicit none
  integer :: ifvcfpout,ifhvccfp,is,  if1011,if3011, ifplane,ngpmx, ngcmx, nbloch,&
       ibas,ic,lxx,nxx,nrx,l,n,k,isx,kdummy, nkdmx,nkqmx,lmax,nkdest,nkrest,ngp,ngc,nlxx,i,lnjcg,lnxcg, &
       nkd,nkq ,ibas1,ibas2,nlx1,nlx2, iqibz,ir,ig1,n1,n2, ngb,nev,nmx,iqx,ipl1,ipl2,igx1,igx2,&
       igc,igc0,ifgb0vec,ifgb0vec1,ix, iy, iqxini, iqxend,imode, ngc0, ifvcfporg,nqbz_in,&
       ifprodmt,nl_r,lx_,nxx_r,nxdim,ibl1,nn,no,ngbnew, nmatch,ifpmatch,nmatch_q,ifpmatch_q,m,ifpomat,nbln,ibln,ngb_in,nnr,igc2,&
       nnmx ,ngcnn, ifvcoud,idummy,ifiwqfac,iqbz,iqbzx,nnn,ixyz,ifq0p,incwfin, &
       nqnumc,ifiqgc , mpi__iini, mpi__iend
  integer,allocatable :: jcg(:),indxcg(:), lx(:),kmx(:),nblocha(:),nr(:),ificrb(:), &
       nx(:,:),ngvecp(:,:),ngvecc(:,:),ngvecci(:,:,:),iqibzx(:)
  integer,allocatable:: ngvecc0(:,:)
  integer,allocatable:: nx_r(:), ibl(:,:,:,:),imatcho(:),imatchn(:),imatcho_q(:),imatchn_q(:)
  real(8) ::  q(3),p(3),voltot, tripl,alat0,epsx, tol,as,tpiba,qb0(3,3),vol0,rdist0,qdist0,radd,qadd, &
       a0,awald,alat1,tol1,r0,q0,awald0,qg(3),   absqg2,aaa,aaa12, eee,eees, q_org(3),&
       absqq,qqx(3), epsmx,aaaa, sss1,sss2,dnorm, qqq(3),QpGcut_Cou,quu(3)
  real(8),allocatable:: prodmt(:,:,:,:),rdmatch(:,:,:,:), rmax(:), cg(:),rprodx(:,:,:,:),dlv(:,:),qlv(:,:),work(:),ngcn(:), &
       rojb(:,:,:),sgbb(:,:,:,:),aa(:),bb(:),rofit(:),phi(:),psi(:),wqfac(:),qbzwww(:,:),rkpr(:,:,:),rkmr(:,:,:),rofi(:,:), eb(:)
  real(8),parameter::pi  = 4d0*datan(1d0), fpi = 4d0*pi
  complex(8):: pval,pslo,phasex, phasep,img=(0d0,1d0), xxx,trwv
  complex(8),allocatable:: geig(:,:),strx(:,:,:,:),sgpb(:,:,:,:),sgpp(:,:,:,:), fouvb(:,:,:,:),fouvp(:,:,:,:),&
       vcoul0(:,:), s(:,:),sd(:,:),rojp(:,:,:) , vcoulnn(:,:), gbvec(:), vcoul_org(:,:),&
       matp(:),matp2(:),ppmt(:,:,:,:),pmat(:,:),pomat(:,:),zzr(:)
  logical :: checkeig, besseltest=.false.,smbb, wvcc, cmdopt2,cmdopt0,debug=.false. !emptyrun,
  character(20) :: xxt,outs=''
  character(3) :: charnum3
  character(10) :: i2char
  character(128):: vcoudfile
  character(1024):: aaaw
  real(8),external::screenfac
  complex(8),allocatable,target :: oo(:,:)
  complex(8),allocatable,target :: vcoul(:,:) !,hh(:,:),zz(:,:)
  complex(8),pointer:: ppovl(:,:),zz(:,:),hh(:,:)
  call MPI__Initialize()
  call gpu_init(comm) 
  call M_lgunit_init()
!  emptyrun=cmdopt0('--emptyrun')
  if( mpi__root) write(6,"(' mode=0,3,202 (0 and 3 give the same results for given bas)' )")
  if(cmdopt2('--job=',outs)) then; read(outs,*) imode
  elseif( mpi__root ) then       ; read(5,*) imode;   endif
  call MPI__Broadcast(imode) !  write(ixcc,"('.mode=',i4.4)")imode
  call MPI__consoleout('hvccfp0.mode'//charnum3(imode))
  call cputm(stdo)
  if(imode==202 )  then; if(ipr) write(6,*)' hvccfp0: imode=',imode
  elseif(imode==0) then
  elseif(imode==3) then
  else; call rx('hvccfp0: now hvccfp0 support just normal mode=0 3 202')
  endif
  call Genallcf_v3(incwfx=0) !use ForX0 for core in GWIN
  call Readhamindex0()
  call Readhamindex()
  call Read_bzdata()
  call readngmx('QGcou',ngcmx)
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  allocate(ngvecc(3,ngcmx))
  call ReadGWinputKeys() !Readin dataset in GWinput
  ReadinBASFP:block
    allocate(lx(nbas),kmx(nbas),nblocha(nbas),nr(nbas),aa(nbas),bb(nbas),ificrb(nbas),rmax(nbas) )
    do ibas = 1,nbas !! Readin BASFP//atom. The product basis functions.
       ic = ibas !
       open(newunit=ificrb(ibas),file='__BASFP'//char( 48+ic/10 )//char( 48+mod(ic,10)))
       read(ificrb(ibas),"(4i6,2d24.16)") lx(ibas), kmx(ibas), nblocha(ibas), nr(ibas),aa(ibas),bb(ibas)
       rmax(ibas) = bb(ibas)*(exp((nr(ibas)-1)*aa(ibas))-1d0)
    enddo
    lxx = maxval(lx)
    allocate( nx(0:lxx,nbas) )
    do ibas = 1,nbas
       read(ificrb(ibas),"(i5)") nx(0:lx(ibas),ibas)
    enddo
    nxx = maxval(nx)
    nrx = maxval(nr)
    allocate( rprodx(nrx,nxx,0:lxx,nbas) )
    do ibas = 1,nbas
       do l = 0, lx(ibas)
          do n = 1, nx(l,ibas)
             read(ificrb(ibas),"(3i5)"   ) k, kdummy,kdummy
             read(ificrb(ibas),"(d23.15)") (rprodx(i,n,l,ibas),i=1,nr(ibas))
          enddo
       enddo
       close(ificrb(ibas))
    enddo
  endblock ReadinBASFP

  besseltest=.False.
  BesselTestBlock: if(besseltest) then
     if(ipr) write(6,*)' *** TEST case ***  rprodx is given by Bessel.'
     iqx  = 2 ! test G, corresponding <q+G|v|q+G> should be exact. e.g. ig1=1 and  ig1=35 for iqx=2 
     igx1 = 1
     igx2 = 35
     if(ipr) write(6,"(' iqx=',i3,' ig1 ig2=',2i3)") iqx,igx1,igx2
     if(ipr) write(6,"(a)") ' <q+G|v|q+G> for the corresponding iqx ig1 ig2 should be exact!'
     if(ipr) write(6,"(a)") ' See fort.196'
     if(ipr) write(6,"(a)") ' Errors will be from the radial function integrals !!!'
     if(ipr) write(6,"(a)") ' You can slso so similar test from hbasfp0.'
     if(ipr) write(6,"(a)") ' See test1 in basnfp0.'
     deallocate(rprodx,nx)
     tpiba=8.d0*datan(1.d0)/alat
     lx = 4
     nr = nr(1)
     aa = aa(1)
     bb = bb(1)
     lxx = maxval(lx)
     allocate( nx(0:lxx,nbas) )
     kmx= 1
     nx = 2
     nxx = maxval(nx)
     nblocha= nxx *(lxx+1)**2
     nrx = maxval(nr)
     allocate(rprodx(nrx,nxx,0:lxx,nbas),rofit(nrx),phi(0:lxx),psi(0:lxx))
     rofit(1) = 0d0
     do ir = 1, nrx;    rofit(ir) = bb(1)*( exp(aa(1)*(ir-1)) - 1d0);   enddo
     do n = 1, nxx
        if(n==1) ig1 = igx1
        if(n==2) ig1 = igx2
        qg(1:3) = tpiba * (qibz(1:3,iqx)+ matmul(qlat, ngvecci(1:3,ig1,iqx)))
        absqg2  = sum(qg(1:3)**2)
        do ir =1,nrx
           call bessl(absqg2*rofit(ir)**2,lxx,phi,psi)
           do ibas=1,nbas
              do l = 0, lx(ibas)
                 rprodx(ir,n,l,ibas) = phi(l)* rofit(ir) **(l +1 )
              enddo
           enddo
        enddo
     enddo
     do ibas=1,nbas
        do l = 0, lx(ibas) !orthogonalized rprodx.
           rprodx(1:nr(ibas),1,l,ibas)= rprodx(1:nr(ibas),1,l,ibas) + rprodx(1:nr(ibas),2,l,ibas)
           n = 1
           call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas),aa(ibas),bb(ibas),nr(ibas), aaa )
           aaa = 1d0/sqrt(aaa)
           rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
           if(nxx==1) cycle
           n1=1
           n2=2
           call gintxx(rprodx(1,n1,l,ibas),rprodx(1,n2,l,ibas),aa(ibas),bb(ibas),nr(ibas), aaa12 )
           rprodx(1:nr(ibas),n2,l,ibas) = rprodx(1:nr(ibas),n2,l,ibas) - aaa12*rprodx(1:nr(ibas),n1,l,ibas)
           n = 2
           call gintxx(rprodx(1,n,l,ibas),rprodx(1,n,l,ibas),aa(ibas),bb(ibas),nr(ibas), aaa )
           aaa = 1d0/sqrt(aaa)
           rprodx(1:nr(ibas),n,l,ibas)= aaa*rprodx(1:nr(ibas),n,l,ibas)
        enddo
     enddo
  endif BesselTestBlock
  
  CGcoefficienets4realharmonics: block
    ! <LM3|lm1 lm2>.  We have one another GGgenerator, clebsh_t, easier indexing but many zeros.
    ! indx= lm1(lm1-1)/2 + lm2 (Here lm1>lm2. Search cg and indxcg in m_hvccfp0_util)
    ! icg1 = indxcg(indx)
    ! icg2 = indxcg(indx+1)-1
    ! do igc = icg1,icg2
    !    lm3 = jcg(icg)   !note lm3= l3**2 + m+l3 +1
    !   <lm3|lm1 lm2> = cg(icg)
    ! enddo
    use m_scg,only: scg,scg_sizechk
    call scg_sizechk(lmax=lxx, lnjcg=lnjcg,lnxcg=lnxcg) ! return lnjcg,lnxcg for allocation
    allocate(cg(lnjcg),jcg(lnjcg),indxcg(lnxcg))
    call scg(lxx,cg,indxcg,jcg)
  endblock CGcoefficienets4realharmonics

  Mesh4Ewald: block
    tol    = 1d-14
    ! tol was 1d-9 before aug2019.
    ! Sakakibara had a problem for slab model of KF for 18 atoms per cell. --> The lowest eigenvalue of coulomb matrix becomes negative.
    !Apr2021 gomi. we had a problem of posititve definiteness of the Coulomb matrix. Probably because of this routine   ! These default values ok?
    awald0 = 2d0   ! See plat_0
    nkdmx  = 8000  !    800->8000 to treat GdN. !31oct2019
    nkqmx  = 8000  !    nkdmx and nqkmx should be determined automatically in future.
    lmax   = 2*lxx ! lxx or lmax=6 ???
    vol0= abs(tripl(plat,plat(1,2),plat(1,3)))
    as   = awald0
    alat1= alat
    tpiba= 8.d0*datan(1.d0)/alat
    call cross(plat(1,2),plat(1,3),qb0(1,1))
    call cross(plat(1,3),plat(1,1),qb0(1,2))
    call cross(plat(1,1),plat(1,2),qb0(1,3))
    qb0(1:3,1:3) = qb0(1:3,1:3)/vol0
    rdist0=vol0**(1.d0/3.d0)
    qdist0=1.d0/rdist0
    radd=.7*rdist0
    qadd=.7*qdist0
    a0=as/rdist0
    awald=a0/alat
    tol1= tol*alat**(lmax+1)
    allocate(dlv(3,nkdmx), qlv(3,nkqmx), work(max0(nkdmx,nkqmx)) )
    call lctoff(a0,vol0,lmax,tol1,r0,q0)
    nkdest =4.18879*(r0+radd)**3/vol0+.5
    nkrest =4.18879*(q0+qadd)**3*vol0+.5
    call lgen(plat,r0+radd,nkd,nkdmx,dlv,work)
    call lgen(qb0, q0+qadd,nkq,nkqmx,qlv,work)
    if(ipr) write(6,"(/' lattc:  as=',f6.3,'   tol=',1p,e9.2,'   lmax=',i2,'   awald=',0p,f7.4,'   v0=',f10.3/' alat1=',f9.5, &
         '   estimates:   nkd',i6,'   nkr',i6)") as,tol,lmax,awald,vol0,alat1,nkdest,nkrest
    if(ipr) write(6,"('  r0=',f9.4,'   rc=',f9.4,'   radd=',f9.4,'   nkd=', i7)") r0,r0*alat,radd,nkd
    if(ipr) write(6,"('  q0=',f9.4,'   qc=',f9.4,'   qadd=',f9.4,'   nkr=', i7)") q0,q0*tpiba,qadd,nkq
    deallocate(work)
  endblock Mesh4Ewald
  
  qindependentRadialIntegrals:block
    eee= merge(0d0,screenfac(),imode==202)
    if(ipr) write(6,"(' Coulomb is exp(sqrt(-eee)*r)/r. eee=',d13.6,d13.6)")eee
    ! For eps_lmf and epsPP_ mode, even small eee=1d-4 affects to dielectric function near q=0 when eps values is large ~>100 or more.
    ! Thus we set eee=0d0 to avoid this.
    ! bessel and hankel for the expansion of exp(-r/r_0)/r. bessel and hankel is renomarized so that its behaves as r^l and r^{-l-1} near r=0.
    !  rkpr means r^l*r for e=0 (r0c =infinity) case
    allocate(rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas),rofi(nrx,nbas))
    do ibas=1,nbas
       call genjh(eee,nr(ibas),aa(ibas),bb(ibas),lx(ibas), nrx,lxx,rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas))
    enddo
    allocate( rojb(nxx, 0:lxx, nbas), sgbb(nxx,  nxx,  0:lxx, nbas))
    do ibas = 1,nbas !
       call mkjb_4(lxx, lx(ibas),nxx, nx(0:lxx,ibas), aa(ibas),bb(ibas), nr(ibas), nrx,rprodx(1,1,0,ibas), &
            rofi(1,ibas),rkpr(1,0,ibas),rkmr(1,0,ibas),  rojb(1,0,ibas),     sgbb(1,1,0,ibas))
    enddo                              !Onsite integrals rojb=<j(e=0)|B> and sgbb=<B|v(onsite)|B>
  endblock qindependentRadialIntegrals
  call cputm(stdo, 'end of qindependentRadialIntegrals')

  nlxx= (lxx+1)**2
  allocate(ngvecc0(3,ngcmx))
  call readqg('QGcou',[0d0,0d0,0d0],  quu,ngc0, ngvecc0) !coulomb matrix for each q = qibz
  deallocate(ngvecc0)
  nbloch    = sum(nblocha)
  iqxend = nqibz + nq0i+nq0iadd
  if(imode==202) then; iqxini= nqibz + 1
  else;                iqxini = 1
  endif
  if(ipr) write(6,*)'iqxini iqxend=',iqxini,iqxend
  if(abs(sum(qibz(:,1)**2))/=0d0) call rx( 'hvccfp0: We assume sum(q**2)==0d0 but not.')
  call MPI__getRange( mpi__iini, mpi__iend, iqxini, iqxend )

  mainforiqx: do 1001 iqx = mpi__iini, mpi__iend !q=(0,0,0) is omitted!    !write(stdo,"('#### do 1001 start iqx=',5i5)")iqx,nqibz
    if(iqx<=nqibz) then; q= qibz(1:3,iqx)
    else;                q= q0i(1:3,iqx-nqibz)
    endif
    if(imode==202 .AND. abs(sum(q))<1d-8) cycle
    call readqg('QGcou',q,  quu,ngc, ngvecc ) !Get q+G vector
    ngb = nbloch + ngc  
    write(aaaw,'(" do 1001: iq nqibz q ngc =",2i5,3f10.4,i5)') iqx,nqibz,q,ngc
    call cputm(stdo,aaaw)
    allocate( strx(nlxx,nbas,nlxx,nbas), source = (0d0,0d0)) !! strxq: structure factor.
    do ibas1 =1,nbas
      do ibas2 =1,nbas
        p = bas(:,ibas2)-bas(:,ibas1)
        phasep =exp(img*2*pi*sum(q*p))
        nlx1 = (lx(ibas1)+1)**2
        nlx2 = (lx(ibas2)+1)**2
        allocate( s(nlx1,nlx2),sd(nlx1,nlx2)) !kino add sd----but sd is dummy
        call strxq(1,eee,q,p,nlx1,nlx2,nlx1,alat,voltot,awald,nkd,nkq,dlv,qlv,cg,indxcg,jcg, s,sd)
        strx(1:nlx1,ibas1,1:nlx2,ibas2) = fpi*s      !!! *phasep
        deallocate( s,sd )
      enddo
    enddo

    call cputm(stdo,'goto mkjp_4')
    qdependentRadialIntegrals:block
      allocate( rojp(ngc, nlxx, nbas), sgpb(ngc, nxx, nlxx, nbas), fouvb(ngc, nxx, nlxx, nbas))
      do ibas = 1,nbas 
        call mkjp_4(q,ngc, ngvecc, alat, qlat, lxx, lx(ibas),nxx, nx(0:lxx,ibas), bas(1,ibas),aa(ibas),bb(ibas),rmax(ibas), &
             nr(ibas), nrx, rprodx(1,1,0,ibas), eee, rofi(1,ibas), rkpr(1,0,ibas), rkmr(1,0,ibas), &
             rojp(1,1,ibas),               sgpb(1,1,1,ibas),                   fouvb(1,1,1,ibas)) 
      enddo! rojp=<j(e=0)_L|exp(i(q+G)r)>, sgpb=<exp(i(q+G)r)|v(onsite)|B_nL>, fouvb=<exp(i(q+G)r)|v|B_nL>
    endblock qdependentRadialIntegrals

    write(aaaw,ftox)'start vcoulq_4 for iqx=',iqx,' q=',ftof(q(1:3),4),' ngc=',ngc,' nbas=',nbas,' lx=',lxx,' nx=',nxx,'procid=',mpi__rank
    call cputm(stdo,aaaw)
    
    allocate( vcoul(ngb,ngb), source=(0d0,0d0) )
    call vcoulq_4(q, nbloch, ngc, nbas, lx,lxx, nx,nxx, alat, qlat, voltot, ngvecc, strx, rojp,rojb, sgbb,sgpb, fouvb, ngb, &
         bas,rmax, eee, aa,bb,nr,nrx,rkpr,rkmr,rofi, &
         vcoul) !the Coulomb matrix
    deallocate( strx, rojp,sgpb,fouvb)
    
    write(aaaw,ftox)'end of vcoulq_4 for iqx=',iqx,' q=',q(1:3),' ngc=',ngc,' nbas=',nbas,' lx=',lxx,' nx=',nxx,'procid=',mpi__rank
    call cputm(stdo,aaaw)
    if(debug) then
       write(stdo,'(" vcoul trwi=",i6,2d22.14)') iqx,sum([(vcoul(i,i),i=1,nbloch)])
       write(stdo,'("### sum vcoul(1:ngb,      1:ngb) ",2d22.14,2x,d22.14)') sum(vcoul), sum(abs(vcoul))
       write(stdo,'("### sum vcoul(1:nbloch,1:nbloch) ",2d22.14,2x,d22.14)') sum(vcoul(1:nbloch,1:nbloch)),sum(abs(vcoul(1:nbloch,1:nbloch)))
       write(stdo,*)
    endif
    allocate( oo(ngb,ngb)  ,source=(0d0,0d0))
    forall(ipl1=1:nbloch) oo(ipl1,ipl1) = 1d0
    ppovl => oo(nbloch+1:ngb,nbloch+1:ngb) !ppovl is a part of oo
    call mkppovl2(alat,plat,qlat, ngc,  ngvecc, ngc,  ngvecc, nbas, rmax, bas, ppovl)
    call cputm(stdo,' end of mkppovl2')

    Diagonalize_Coulomb_matrix: block
#ifdef __GPU
      use m_lapack, only: zhgv => zhgv_d
#else
      use m_lapack, only: zhgv => zhgv_h
#endif
      integer :: istat
      allocate(eb(ngb))
      nev = ngb
      vcoul=-vcoul ! trick to get decendent order of eigenvalues by zhgv.
      !$acc data copyin(vcoul) copy(oo) copyout(eb)
      istat = zhgv(vcoul, oo, ngb, eb)
      !$acc end data
      deallocate(oo)
      zz => vcoul
    endblock Diagonalize_Coulomb_matrix
    Chkwriteeb: do ipl1=1,ngb
      if(ipl1==1 .and.ipr) write(6,*)' --- goto eigen check1 --- '
      if(ipl1==11.and.ipr) write(6,*)' ... '
      if(ipl1>10 .AND. ipl1<nev-5) cycle
      if(ipr) write(6,'(i4,d23.16)')ipl1,-eb(ipl1)
    enddo Chkwriteeb
    if(ipr) write(6,"(' nev ngv q=',2i5,3f10.6)")nev,ngb,q
    !! -eb should be positive definite. However, we have one (or a few?) negative ones.
    !! TK think no problem to set eb=0 when -eb is negative.
    !! But this is a temporaly fix or better manner to calcuate coulomb matrix. strxq may be needed to be replaced
    !do i=1,nev
    !  if(eb(i)>0 .AND. keeppositivecou) then
    !     if(ipr) write(6,"(a,d13.5)")'KeepPositiveCou enforce : -eb<0 --> eb=0',eb(i)
    !     eb(i)=0d0
    !  endif
    !enddo
    eb = merge(0d0,eb,-eb<0 .AND. keeppositivecou) !'KeepPositiveCou enforce : -eb<0 --> eb=0',eb(i)
    if(debug) then !we expect vcoul is reserved. 
      if(ipr) write(6,*)
      if(ipr) write(6,'(" eig0 must be equal to the largest =", 2d24.16)') sum(  dconjg(zz(1:ngb,1))*matmul( vcoul(1:ngb,1:ngb),zz(1:ngb,1)))
 !     if(ipr) write(6,'(" zz norm check=",d24.16)')    sum( dconjg(zz(1:ngb,1))*matmul(oox,zz(1:ngb,1)) )
      if(ipr) write(6,*)
      if(ipr) write(6,'(" --- vcoul(exact)=",d14.6," absq2=",d24.16)') fpi*voltot/(sum(tpiba**2*q(1:3)**2)-eee), (sum(tpiba**2*q(1:3)**2)-eee)
      if(ipr) write(6,'(" --- vcoul(cal ) =",2d14.6)') sum( dconjg(zz(1:ngb,1))*matmul( vcoul(1:ngb,1:ngb),zz(1:ngb,1)) )*voltot
     endif 
    open(newunit=ifvcoud,file=trim('__Vcoud.'//i2char(iqx)),form='unformatted') !  !! Vcoud file, which contains E(\nu,I), given in qibzPRB81,125102
    write(ifvcoud) ngb
    write(ifvcoud) q
    write(ifvcoud) -eb
    write(ifvcoud) zz !=Enu
    close(ifvcoud)
    deallocate(eb,vcoul)
    call cputm(stdo,' end of do 1001')
1001 enddo mainforiqx
  deallocate(ngvecc)
  call cputm(stdo,'end of hvccfp0')
  if(imode==202) call rx0( ' OK! hvccfp0 imode=202 only for Q0P')
  if(imode==0)   call rx0( ' OK! hvccfp0 imode=0')
  if(imode==3)   call rx0( ' OK! hvccfp0 imode=3')
end subroutine

subroutine MPI__getRange( mpi__indexi, mpi__indexe, indexi, indexe )
  use m_mpi,only: mpi__size,mpi__rank
  implicit none
  integer, intent(out) :: mpi__indexi, mpi__indexe
  integer, intent(in)  :: indexi, indexe
  integer, allocatable :: mpi__total(:)
  integer              :: total
  integer :: p
  allocate( mpi__total(0:mpi__size-1) )
  total = indexe-indexi+1
  mpi__total(:) = total/mpi__size
  do p=1, mod(total,mpi__size)
     mpi__total(p-1) = mpi__total(p-1) + 1
  end do
  mpi__indexe=indexi-1
  do p=0, mpi__rank
     mpi__indexi = mpi__indexe+1
     mpi__indexe = mpi__indexi+mpi__total(p)-1
  end do
  deallocate(mpi__total)
end subroutine MPI__getRange
end module
